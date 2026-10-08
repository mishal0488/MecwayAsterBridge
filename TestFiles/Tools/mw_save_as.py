# -*- coding: utf-8 -*-
# Copyright (C) 2026 trednuK
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Mecway-script (IronPython 2.7): het geopende model opslaan zoals Mecway
# het zelf schrijft (eigen XML-opmaak, volgorde, eenheden).
# gebruik: mecway.exe <in.liml> script mw_save_as.py
#   (Mecway-handleiding par. 16.1: "script opens the file <filename>, runs
#   the script <script filename>, then exits"; geen solve, dus ook voor
#   analyses die Mecway's eigen solvers niet kunnen rekenen.)
# Doel: omgevingsvariabele MAB_MW_OUT, anders <in>_mw.liml naast het model.
# Gebruikt door TestFiles/Verificatie/run_verificatie.py --mecway-opslaan.
import os

src = mw.file_name()
dst = os.environ.get("MAB_MW_OUT") or (os.path.splitext(src)[0] + "_mw.liml")
mw.save_as(dst)
