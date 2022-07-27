# Docker for pyopenms-qc
- `docker-compose up`

## How to run
- `cd /Users/gongm/Documents/projects/Denver_HEU_Pilot/Denver_HEU_07112022/docker/pyopenms-qc/`
- Build: `docker-compose build`
- Start: `docker-compose up`
- Stop: `docker-compose down`


## Port
- The port is changed from '8888:8888' to '8880:8888'. So when you see the below message 

```
openms-qc  |     Or copy and paste one of these URLs:
openms-qc  |         http://dcabd27899f5:8888/?token=f04a4f8d6f8e79334065ed9b8afc316de5e289fbbe05397d
openms-qc  |      or http://127.0.0.1:8888/?token=f04a4f8d6f8e79334065ed9b8afc316de5e289fbbe05397d
```

- You need to change the port from `8888` to `8880` in the link.

## Fix
- Somehow the openms cannot install in the newest jupyter/scipy
- So, I use the older version: `jupyter/scipy-notebook:2021-11-10`