FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive MOJO_MODE=production MOJO_REVERSE_PROXY=1
RUN apt-get update && apt-get install -y --no-install-recommends \
      perl libbio-perl-perl libmojolicious-perl libdbi-perl libdbd-sqlite3-perl sqlite3 \
      python3 python3-numpy python3-pandas python3-sklearn python3-joblib \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY . /app
# Retrain the model and rebuild the database inside the image so they match its library versions.
RUN python3 student3/train_extratrees.py data/final/training_balanced_props.tsv --outdir data/final/model \
 && python3 student3/predict_et.py data/final/model/model_et.joblib data/final/discovery_marine_props.tsv > data/final/discovery_marine_scores_et.tsv \
 && perl student4/build_db.pl data/final > data/final/marine.sql \
 && sqlite3 data/final/marine.db < data/final/marine.sql
EXPOSE 10000
CMD ["sh", "-c", "perl student4/app.pl daemon -l http://0.0.0.0:${PORT:-10000}"]
