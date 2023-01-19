#!/usr/bin/python2

import sys
import os
import subprocess
from subprocess import getstatusoutput
from subprocess import getoutput
import datetime
import math

from parser_utils import *

def getProxy():
    stat,out = subprocess.getstatusoutput("voms-proxy-info -e --valid 5:00")
    if stat:
        raise Exception("voms proxy not found or validity  less than 5 hours:\n%s" % out)
    stat,out = subprocess.getstatusoutput("voms-proxy-info -p")
    out = out.strip().split("\n")[-1] # remove spurious java info at ic.ac.uk
    if stat:
        raise Exception("Unable to voms proxy:\n%s" % out)
    proxy = out.strip("\n")
    return proxy


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def lxbatchSubmitJob (run, firstspill, nfiles, path, cfg, outdir, queue, job_dir, dryrun):
    lastspill = firstspill+nfiles-1
    jobname = job_dir+'/H4Reco_'+queue+'_'+run+'_{fs}-{ls}.sh'.format(fs=firstspill, ls=lastspill)
    jobtar = job_dir+'/job.tar'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('cp '+jobtar+' ./ \n')
    f.write ('tar -xf job.tar \n')
    f.write ('source scripts/setup.sh \n')
    f.write ('make -j 2 \n')
    f.write ('cp '+path+'/'+cfg+' job.cfg \n\n')
    f.write ('bin/H4Reco job.cfg {r} {fs} {nf}\n\n'.format(r=run, fs=firstspill, nf=nfiles))
    f.write ('cp ntuples/*run{r}_spills{fs}-{ls}.root {o}\n'.format(r=run, fs=firstspill, ls=lastspill, o=outdir))
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+job_dir+'; bsub -q ' + queue + ' ' + '-u ' + os.environ['USER'] + '@cern.ch ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def htcondorSubmitJob(runs, path, cfg, outdir, queue, job_dir, dryrun, notar, version, outname):
    jobname = job_dir+'/H4Reco_condor'
    jobtar = job_dir+'/job.tar'
    #---H4Reco script
    fsh = open (jobname+'.sh', 'w')
    fsh.write ('#!/bin/sh' + '\n\n')
    fsh.write ('declare -a runs=('+' '.join([str(run) for run in runs])+')\n\n')
    fsh.write ('cp '+jobtar+' ./ \n')
    if not notar:
        fsh.write ('tar -xf job.tar \n')
        fsh.write ('source scripts/setup.sh \n')
        fsh.write ('make -j 2 \n')
        fsh.write ('cp '+path+'/'+cfg+' job.cfg \n\n')
        fsh.write ('bin/H4Reco job.cfg ${runs[${1}]}\n\n')
    else:
        fsh.write ('cd '+path+' \n')
        fsh.write ('source scripts/setup.sh \n')
        fsh.write ('bin/H4Reco '+cfg+' ${runs[${1}]}\n\n')
    fsh.write ('cp '+outname+'/*${runs[${1}]}.root '+outdir+'\n')
    fsh.close ()
    #---HTCondor submit file
    fsub = open (jobname+'.sub', 'w')
    fsub.write('+JobFlavour = "'+queue+'"\n\n')
    fsub.write('executable  = '+jobname+'.sh\n')
    fsub.write('arguments   = $(ProcId)\n')
    fsub.write('output      = '+job_dir+'/output/h4reco.$(ClusterId).$(ProcId).out\n')
    fsub.write('error       = '+job_dir+'/output/h4reco.$(ClusterId).$(ProcId).err\n')
    fsub.write('log         = '+job_dir+'/log/h4reco.$(ClusterId).log\n\n')
    fsub.write('x509userproxy = '+getProxy()+' \n\n')
    fsub.write('max_retries = 3\n')
    fsub.write('queue '+str(len(runs))+'\n')
    fsub.close()
    #---submit
    getstatusoutput('chmod 755 '+jobname+ '*')
    getstatusoutput('mkdir -p '+job_dir+'/output')
    getstatusoutput('mkdir -p '+job_dir+'/log')
    if not dryrun:
        ret = getstatusoutput('condor_submit '+jobname+'.sub')
        print(ret)

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def htcondorSubmitJobs (run, nfiles, path, cfg, outdir, queue, job_dir, dryrun, notar,version, outname):
    # prepare submission
    nspills = getNumberOfSpills(run, path, cfg)
    print('run {r} with {s} spills'.format(r=run, s=nspills))
    if nfiles <= nspills:
        njobs = int(math.ceil(nspills / float(nfiles)))
    else:
        njobs = 1
        nfiles = nspills
    print('submitting {j} jobs for run {r}'.format(j=njobs, r=run))

    jobname = job_dir+'/H4Reco_condor_run{r}'.format(r=run)
    jobtar = job_dir+'/job.tar'
    #---H4Reco script
    fsh = open (jobname+'.sh', 'w')
    fsh.write ('#!/bin/sh' + '\n\n')
    if not notar:
        fsh.write ('cp '+jobtar+' ./ \n')
        fsh.write ('tar -xf job.tar \n')
        fsh.write ('source scripts/setup.sh \n')
        fsh.write ('make -j 2 \n')
    else: 
        fsh.write ('cd '+path+' \n')
        fsh.write ('source scripts/setup.sh \n')

    fsh.write ('firstspill=$((${1}*'+str(nfiles)+'+1))\n')
    fsh.write ('lastspill=$(((${1}+1)*'+str(nfiles)+'))\n\n')
    fsh.write ('bin/H4Reco '+cfg+' {r} $firstspill {nf}\n\n'.format(r=run, nf=nfiles))
    fsh.write ('mv '+outname+'/{r}/$firstspill.root {o}\n'.format(v=version, r=run, o=outdir))
    fsh.close ()
    #---HTCondor submit file
    fsub = open (jobname+'.sub', 'w')
    fsub.write('+JobFlavour = "'+queue+'"\n\n')
    fsub.write('executable  = '+jobname+'.sh\n')
    fsub.write('arguments   = $(ProcId)\n')
    fsub.write('output      = '+job_dir+'/output/h4reco.$(ClusterId).$(ProcId).out\n')
    fsub.write('error       = '+job_dir+'/output/h4reco.$(ClusterId).$(ProcId).err\n')
    fsub.write('log         = '+job_dir+'/log/h4reco.$(ClusterId).log\n\n')
    fsub.write('max_retries = 3\n')
    fsub.write('queue {}\n'.format(njobs))
    fsub.close()
    #---submit
    getstatusoutput('chmod 755 '+jobname+ '*')
    getstatusoutput('mkdir -p '+job_dir+'/output')
    getstatusoutput('mkdir -p '+job_dir+'/log')
    if not dryrun:
        ret = getstatusoutput('condor_submit '+jobname+'.sub')
        print(ret)

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def herculesSubmitJob (run, firstspill, nfiles, path, cfg, outdir, queue, job_dir, dryrun):
    lastspill = firstspill+nfiles-1
    jobname = job_dir+'/H4Reco_'+queue+'_'+run+'_{fs}-{ls}.sh'.format(fs=firstspill, ls=lastspill)
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('mkdir -p /gwpool/users/$USER/pool/'+run+'/ \n')
    f.write ('cd /gwpool/users/$USER/pool/'+run+'/ \n')
    f.write ('wget https://github.com/simonepigazzini/H4Analysis/archive/master.zip \n')
    f.write ('unzip master.zip \n')
    f.write ('cd H4Analysis-master/ \n\n')
    f.write ('cp '+path+'/ntuples/Template*.root ./ntuples/ \n')
    f.write ('cp '+path+cfg+' job.cfg \n')
    f.write ('source scripts/setup.sh \n')
    f.write ('make -j 2\n\n')
    f.write ('bin/H4Reco job.cfg {r} {fs} {nf}\n\n'.format(r=run, fs=firstspill, nf=nfiles))
    f.write ('cp ntuples/*run{r}_spills{fs}-{ls}.root {o}\n'.format(r=run, fs=firstspill, ls=lastspill, o=outdir))
    f.write ('cd /gwpool/users/$USER/pool/ \n')
    f.write ('rm -r '+run+' \n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+job_dir+'; qsub -q ' + queue + ' ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
def getPath2data (path, cfg, whichpath = 'path2data'):
    # parse cfg file to find the path parameter
    path2data = ''
    with open(path+'/'+cfg) as cfgfile:
        for line in cfgfile:
            print (line)
            words = line.split()
            if len(words) == 2:
                if whichpath in words[0]:
                    path2data = words[1]
                    break
    return path2data 

def getNumberOfSpills (run, path, cfg):
    nspills = 0
    # parse cfg file to find the path2data parameter
    path2data = getPath2data(path, cfg)
    if not len(path2data):
        with open(path+'/'+cfg) as cfgfile:
            for line in cfgfile:
                words = line.split()
                if len(words) == 2:
                    if 'importCfg' in words[0]:
                        path2data = getPath2data(path, words[1])
                        break
                        
    datafiles = [name for name in os.listdir(path2data+'/'+run+'/') if os.path.isfile(path2data+'/'+run+'/'+name)]

    ## find the number of files in the directory and take this as the number of spills
    #nspills = len(datafiles)

    # use the spill number of the last file as the number of spills assuming a file name like "spill number".root
    spills = [int(datafile.replace('.root', '')) for datafile in datafiles]
    spills.sort()
    nspills = spills[-1]
    return nspills
  
def resubmitFailed(runs, outdir):
    
    successful_runs = getoutput("ls -lhrt "+outdir+"| sed 's:^.*_::g' | sed 's:.root::g'")
    
    successful_runs = set(successful_runs.split('\n'))
    
    return [run for run in runs if run not in successful_runs]    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'submit H4Reco step to lsf')
    parser.add_argument('-r', '--runs' , action=customAction, help='run to be processed, either list or file')
    parser.add_argument('--spills-per-job' , dest='spillsperjob', default = -1, type=int, help='number of spills per job')
    parser.add_argument('-q', '--queue' , default = 'longlunch', help='batch queue/condor flavour (def: longlunch)')
    parser.add_argument('-s', '--storage' , default = '/store/group/dpg_ecal/alca_ecalcalib/ECALTB_H4_Fall2015/', help='storage path')
    parser.add_argument('-v', '--version' , default = 'v1', help='production version')
    parser.add_argument('-c', '--cfg' , default = '../cfg/H4DAQ_base.cfg', help='production version')
    parser.add_argument('--dryrun' , action="store_true", default=False, help='do not submit the jobs, just create them')
    parser.add_argument('--batch' , default='condor', help='batch system to use')
    parser.add_argument('--resub' , action="store_true", default=False, help='resubmit failed jobs')
    parser.add_argument('--notar' , action="store_true", default=False, help='do not create tarball')
    
    args = parser.parse_args ()

    ## check ntuple version
    stageOutDir = args.storage+'/ntuples_'+args.version+'/'
    stageOutDir = stageOutDir.replace('//', '/')
    
    if args.batch == 'lxbatch':
        if getoutput('ls '+stageOutDir) == "":
            print("ntuples version "+args.version+" directory on eos already exist! no jobs created.")
            exit(0)
    getstatusoutput('mkdir -p '+stageOutDir)    
    
    ## job setup
    local_path = getoutput('pwd')
    date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    job_dir = local_path+"/"+date+"_ntuples_"+args.version
    print ("Jobs directory: " + job_dir)
    getstatusoutput('mkdir -p '+job_dir)
    if local_path.endswith('scripts'):
        local_path = local_path[:-len('scripts')]

    if args.cfg.startswith('../'):
        args.cfg = args.cfg[len('../'):]
    
    if len(args.runs) == 1 and os.path.isfile(args.runs[0]):
        runs_file = open(args.runs[0], 'r')
        args.runs  = []
        if runs_file:
            for run in runs_file:
                args.runs.append(run.rstrip())
    outname = getPath2data(local_path, args.cfg, 'outNameSuffix')
    if not len(outname):
        with open(local_path+'/'+args.cfg) as cfgfile:
            for line in cfgfile:
                words = line.split()
                if len(words) == 2:
                    if 'importCfg' in words[0]:
                        path2data = getPath2data(local_path, words[1], 'outNameSuffix')
                        break
    
    ## resubmit failed
    if args.resub:
        args.runs = resubmitFailed(args.runs, stageOutDir)

    ## create jobs
    if not args.notar: getstatusoutput('tar --exclude-vcs --exclude="20*_ntuples*" -cjf '+job_dir+'/job.tar -C '+local_path+' .')
    print('submitting', len(args.runs), 'jobs to queue', args.queue)

    if args.batch == 'condor':
        if args.spillsperjob > 0:
            for run in args.runs:
                htcondorSubmitJobs(run, args.spillsperjob, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun, args.notar, args.version, outname)
        else:
            htcondorSubmitJob(args.runs, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun, args.notar, args.version, outname)
    else:
        for run in args.runs:
            firstspill = 1
            nspills = getNumberOfSpills(run, local_path, args.cfg)
            print('requested to submit run {r} containing {s} spills'.format(r=run, s=nspills))
            nfiles = args.spillsperjob
            if nfiles < 1 or  nfiles > nspills:
                nfiles = nspills

            jobctr = 0
            while firstspill <= nspills:
                if args.batch == 'lxbatch':
                    lxbatchSubmitJob(run, firstspill, nfiles, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun)
                elif args.batch == 'hercules':
                    herculesSubmitJob(run, firstspill, nfiles, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun)
                firstspill += nfiles
                jobctr += 1
            if not args.dryrun:
                print('submitted {j} jobs to {b}'.format(j=jobctr, b=args.batch))
