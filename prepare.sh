#! /bin/bash

./removecomments.sed $1 | sed -e "/^#line/d" -e "/^$/d" -e "/^;/d" | perl -0pi -e 's/\\\n#/\n#/sg' | indent -kr -i8 -ts8 -sob -l80 -ss -bs > temp.c;
mv temp.c $1;
