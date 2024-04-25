This container hosts a built version of HC.

docker run -it --rm -v $HOME/hc:/home/hc_user/work geodynamics/hc

This command will start the HC docker image and give you terminal access. Any changes made in the /home/hc_user/work directory will be reflected on the host machine at home/hc.
