# Minimal Docker image for MultiVirusConsensus using Alpine base
FROM alpine:latest

# install MultiVirusConsensus and dependencies
RUN apk update && \
    apk add --no-cache autoconf automake bash bzip2-dev curl-dev g++ gcc make musl-dev python3 xz-dev zlib-dev && \

    # install htslib (needed for ViralConsensus)
    wget -qO- "https://github.com/samtools/htslib/releases/download/1.24/htslib-1.24.tar.bz2" | tar -xj && \
    cd htslib-* && \
    ./configure && \
    make install && \
    cd .. && \
    rm -rf htslib-* && \
    ldconfig && \

    # install samtools
    wget -qO- "https://github.com/samtools/samtools/releases/download/1.24/samtools-1.24.tar.bz2" | tar -xj && \
    cd samtools-* && \
    ./configure --without-curses && \
    make install && \
    cd .. && \
    rm -rf samtools-* && \

    # install minimap2
    wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.31.tar.gz" | tar -zx && \
    cd minimap2-* && \
    make && \
    mv minimap2 /usr/local/bin/ && \
    cd .. && \
    rm -rf minimap2-* && \

    # install ViralConsensus
    wget -qO- "https://github.com/niemasd/ViralConsensus/archive/refs/tags/1.0.4.tar.gz" | tar -zx && \
    cd ViralConsensus-* && \
    make && \
    mv viral_consensus /usr/local/bin/ && \
    cd .. && \
    rm -rf ViralConsensus-* && \

    # install MultiVirusConsensus
    wget -q -O /usr/local/bin/MultiVirusConsensus.py "https://github.com/niemasd/MultiVirusConsensus/raw/refs/heads/main/MultiVirusConsensus.py" && \
    chmod a+x /usr/local/bin/MultiVirusConsensus.py && \

    # clean up
    rm -rf /tmp/*
