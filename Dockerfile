FROM python:3.10.4-slim-buster
RUN pip install --upgrade pip
RUN pip install --upgrade numpy
RUN pip install --upgrade requests

ADD . /

ENTRYPOINT ["python", "irisplex.py"]
CMD [ "-V" ]
