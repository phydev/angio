FROM ubuntu:20.04

MAINTAINER João Braz Simões <joaosbraz@gmail.com>

RUN mkdir /code

RUN apt-get update

RUN apt-get -y install gfortran liblapack3 libatlas3-base git

WORKDIR /code

VOLUME ["/code"]