FROM ubuntu:18.04 AS s3quilt

RUN apt-get update && apt-get install -y curl python3-pip

ARG GO_VERSION=1.19.2
RUN curl -Lo go${GO_VERSION}.linux-amd64.tar.gz https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz
RUN tar -C . -xzf go${GO_VERSION}.linux-amd64.tar.gz && rm go${GO_VERSION}.linux-amd64.tar.gz
ENV GOPATH /go
ENV PATH $GOPATH/bin:$PATH
RUN /go/bin/go version
WORKDIR /go/src/s3quilt
COPY go.mod go.sum ./

RUN go mod download
RUN go install golang.org/x/tools/cmd/goimports
RUN go install github.com/go-python/gopy

RUN pip3 install pybindgen setuptools wheel

COPY s3quilt.go ./
RUN gopy build -vm=python3 --output s3quilt .
COPY __init__.py s3quilt/
RUN cp -r s3quilt/ /usr/lib/python3/dist-packages/s3quilt/

WORKDIR /root

FROM ubuntu:18.04
RUN apt-get update && apt-get install  -y python3-pip
COPY --from=s3quilt /usr/lib/python3/dist-packages/s3quilt/ /usr/lib/python3/dist-packages/s3quilt/