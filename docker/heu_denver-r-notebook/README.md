# Docker for exp1-notebook
- `docker-compose up`

## How to run
- locate your docker folder
- Build: `docker-compose build`
- Start: `docker-compose up`
- Stop: `docker-compose down`


## Docker images
- The `Dockerfile` is the current docker image for running the analyses


## docker-compose YAML file


### Appendix
- Kill all processes `docker container kill $(docker container ls -q)` 
- Remove images `docker rmi -f 3c28caa55afe`
- BUild docker from Dockerfile `docker build <path2dockerFileFolder>`
  - retag your docker image `tag 059b73abd843 gongm/exp1-r-jupyter:latest`
- To push your docker image to docker hub `docker push gongm/exp1-r-jupyter:latest`
- If you have a file name of the yaml file which is not as default `docker-compose.yaml`
  - Start: `docker-compose -f ./exp1-r-docker-compose.yaml up`
  - Stop: `docker-compose -f ./exp1-r-docker-compose.yaml down` 
  - It is somehow difficult to use this with `docker-compose build` so I still use default name in the current practice.