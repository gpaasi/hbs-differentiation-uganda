# Docker

Build:
  docker build -t ugs-hbs-malaria-ibd .

Run:
  docker run --rm -it -v $(pwd):/project ugs-hbs-malaria-ibd
