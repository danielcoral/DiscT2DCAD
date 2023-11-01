#!/bin/bash
cd /ludc/Home/daniel_c/
source .bashrc
mamba activate depict_env
projects/DVA/Tools/DEPICT/src/python/depict.py projects/DVA/DiscordantT2Dcomp/scripts/08b_Concordant.cfg
projects/DVA/Tools/DEPICT/src/python/depict.py projects/DVA/DiscordantT2Dcomp/scripts/08b_Discordant.cfg