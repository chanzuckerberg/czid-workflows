package s3quilt

import (
	"context"
	"fmt"
	"sync"
	"sync/atomic"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"

	"github.com/aws/aws-sdk-go-v2/service/s3"
)

type chunk struct {
	s3Start uint64
	length  uint64
	idx     int
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
}

func (d *download) downloadChunk(b chunk) error {
	buffer := make([]byte, b.length)
	rangeHeader := fmt.Sprintf("%d-%d", b.s3Start, b.s3Start+b.length)
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
	d.result[b.idx] = string(buffer)
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
		concurrency: 100,
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
