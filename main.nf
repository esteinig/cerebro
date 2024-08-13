#!/usr/bin/env nextflow

/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Cerebro: metagenomic and -transcriptomic diagnostics for clinical production environments
Copyright (C) 2023  Eike Steinig

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

params.team = "CNS"
params.stage_interval = "5s"

process StageFilesCerebroFS {

    input:
    val interval

    output:
    tuple env("TEST"), path("*.txt"), optional: true

    """
    cerebro-fs --team $params.team list > list.txt
    TEST=\$(head -1 list.txt)
    """
}

workflow production {

    Channel.interval(params.stage_interval) | CerebroStage | view

}