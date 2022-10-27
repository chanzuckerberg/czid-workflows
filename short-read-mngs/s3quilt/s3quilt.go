package s3quilt

import (
	"context"
	"fmt"
	"io"
	"os"
	"sync"
	"sync/atomic"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"

	"github.com/aws/aws-sdk-go-v2/service/s3"
)

type chunk struct {
	s3Start    uint64
	length     uint64
	idx        int
	localStart int64
}

type download struct {
	concurrency int
	s3Client    *s3.Client
	wg          *sync.WaitGroup
	err         error
	errored     uint32
	bucket      *string
	key         *string
	result      []string
	outputFile  *os.File
}

type offsetWriter struct {
	offset int64
	file   *os.File
}

func (oW *offsetWriter) Write(p []byte) (int, error) {
	n, err := oW.file.WriteAt(p, oW.offset)
	oW.offset += int64(n)
	return n, err
}

func (d *download) downloadChunk(b chunk) error {
	rangeHeader := fmt.Sprintf("%d-%d", b.s3Start, b.s3Start+b.length)
	resp, err := d.s3Client.GetObject(context.Background(), &s3.GetObjectInput{
		Bucket: d.bucket,
		Key:    d.key,
		Range:  &rangeHeader,
	})
	if err != nil {
		return err
	}

	if d.outputFile != nil {
		writer := offsetWriter{offset: b.localStart, file: d.outputFile}
		_, err := io.Copy(&writer, resp.Body)
		return err
	}

	bytes, err := io.ReadAll(resp.Body)
	if err != nil {
		return err
	}
	d.result[b.idx] = string(bytes)
	return nil
}

func (d *download) downloader(byteRanges <-chan chunk) {
	for byteRange := range byteRanges {
		if d.err == nil {
			if err := d.downloadChunk(byteRange); err != nil {
				// only set the error if it's the first error
				if atomic.CompareAndSwapUint32(&d.errored, 0, 1) {
					d.err = err
				}
			}
		}
		d.wg.Done()
	}
}

func DownloadChunks(bucket string, key string, starts []uint64, lengths []uint64) ([]string, error) {
	cfg, err := config.LoadDefaultConfig(context.Background(), func(lo *config.LoadOptions) error {
		lo.Region = "us-west-2"
		lo.Credentials = aws.AnonymousCredentials{}
		return nil
	})
	if err != nil {
		return []string{}, err
	}

	s3Client := s3.NewFromConfig(cfg)
	chunks := make(chan chunk, len(starts))
	d := download{
		concurrency: 50,
		s3Client:    s3Client,
		bucket:      &bucket,
		key:         &key,
		wg:          &sync.WaitGroup{},
		err:         nil,
		errored:     0,
		result:      make([]string, len(starts)),
	}

	for w := 1; w <= d.concurrency; w++ {
		go d.downloader(chunks)
	}

	for idx, length := range lengths {
		d.wg.Add(1)
		chunks <- chunk{s3Start: starts[idx], length: length, idx: idx}
	}
	close(chunks)
	d.wg.Wait()
	return d.result, d.err
}

func DownloadChunksToFile(bucket string, key string, outputFilePath string, starts []uint64, lengths []uint64) error {
	cfg, err := config.LoadDefaultConfig(context.Background(), func(lo *config.LoadOptions) error {
		lo.Region = "us-west-2"
		lo.Credentials = aws.AnonymousCredentials{}
		return nil
	})
	if err != nil {
		return err
	}

	s3Client := s3.NewFromConfig(cfg)

	outputFile, err := os.Create(outputFilePath)
	if err != nil {
		return err
	}
	defer outputFile.Close()

	chunks := make(chan chunk, len(starts))
	d := download{
		concurrency: 50,
		s3Client:    s3Client,
		bucket:      &bucket,
		key:         &key,
		wg:          &sync.WaitGroup{},
		err:         nil,
		errored:     0,
		outputFile:  outputFile,
	}

	for w := 1; w <= d.concurrency; w++ {
		go d.downloader(chunks)
	}

	var totalBytes int64 = 0
	for idx, length := range lengths {
		d.wg.Add(1)
		chunks <- chunk{s3Start: starts[idx], length: length, idx: idx, localStart: totalBytes}
		totalBytes += int64(length)
	}
	close(chunks)
	d.wg.Wait()
	return d.err
}
