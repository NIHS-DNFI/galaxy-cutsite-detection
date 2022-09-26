# Galaxy for Cut Site Detection
# VERSION 1.0

FROM quay.io/bgruening/galaxy:20.09
LABEL maintainer="Yuh Shiwa, youyuh48@gmail.com"

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_BRAND="Cut Site Detection" \
  GALAXY_CONFIG_CONDA_AUTO_INSTALL=True

COPY config/startup2.sh /usr/bin/startup2.sh
COPY config/my_tool_list.yml $GALAXY_ROOT/my_tool_list.yml
COPY config/my_tool_setup.sh $GALAXY_ROOT/my_tool_setup.sh
COPY config/tool_conf.xml $GALAXY_ROOT/config/tool_conf.xml.sample
COPY config/workflows/* $GALAXY_ROOT/workflows/

COPY siteseq/siteseq.xml /local_tools/siteseq.xml
COPY siteseq/siteseq.py /local_tools/siteseq.py

RUN apt update && apt install -y netcat-openbsd

RUN chmod +x /usr/bin/startup2.sh && \
  chmod +x $GALAXY_ROOT/my_tool_setup.sh

# Install local tools
RUN install-tools $GALAXY_ROOT/my_tool_list.yml && \
  /tool_deps/_conda/bin/conda clean --tarballs --yes && \
  rm -rf /tool_deps/_conda/pkgs

# Install workflow
ENV GALAXY_CONFIG_TOOL_PATH=/galaxy-central/tools/
ENV PATH="${PATH}:/tool_deps/_conda/bin"
RUN startup_lite && \
    galaxy-wait && \
    workflow-install --publish --workflow_path $GALAXY_ROOT/workflows/ \
    -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_EMAIL -p $GALAXY_DEFAULT_ADMIN_PASSWORD

# Prevent changing owership of conda directory by root
USER galaxy
RUN $GALAXY_ROOT/my_tool_setup.sh
USER root

# Prevent deleting installed tool files by docker's volume function
RUN cp -r /galaxy-central/database/shed_tools /local_tools

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup2.sh"]
