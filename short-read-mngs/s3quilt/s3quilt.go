package s3quilt

import (
	"context"
	"fmt"
	"log"
	"os"
	"sync"

	"github.com/aws/aws-sdk-go-v2/config"

	"github.com/aws/aws-sdk-go-v2/service/s3"
)

type chunk struct {
	s3Start    uint64
	localStart uint64
	length     uint64
}

type download struct {
	concurrency int
	s3Client    *s3.Client
	wg          *sync.WaitGroup
	outputFile  *os.File
	bucket      *string
	key         *string
}

func (d *download) downloadChunk(b chunk) error {
	buffer := make([]byte, b.length)
	rangeHeader := fmt.Sprintf("%d-%d", b.s3Start, b.s3Start+b.length-1)
	resp, err := d.s3Client.GetObject(context.Background(), &s3.GetObjectInput{
		Bucket: d.bucket,
		Key:    d.key,
		Range:  &rangeHeader,
	})
	if err != nil {
		return err
	}
	_, err = resp.Body.Read(buffer)
	if err != nil {
		return err
	}

	_, err = d.outputFile.WriteAt(buffer, int64(b.localStart))
	return err
}

func (d *download) downloader(byteRanges <-chan chunk) {
	for byteRange := range byteRanges {
		// TODO: handle error
		err := d.downloadChunk(byteRange)
		if err != nil {
			log.Fatal(err.Error())
		}
		d.wg.Done()
	}
}

func DownloadChunks(region string, bucket string, key string, outputFilePath string, starts []uint64, lengths []uint64) error {
	cfg, err := config.LoadDefaultConfig(context.Background(), func(lo *config.LoadOptions) error {
		lo.Region = region
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
		concurrency: 100,
		s3Client:    s3Client,
		bucket:      &bucket,
		key:         &key,
		outputFile:  outputFile,
		wg:          &sync.WaitGroup{},
	}

	for w := 1; w <= d.concurrency; w++ {
		go d.downloader(chunks)
	}

	var i uint64 = 0
	for idx, length := range lengths {
		d.wg.Add(1)
		chunks <- chunk{s3Start: starts[idx], length: length, localStart: i}
		i += length
	}
	close(chunks)
	d.wg.Wait()
	return nil
}
