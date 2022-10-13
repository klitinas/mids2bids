#!/usr/bin/env python

# Functions toward converting UMICH ('MIDS') MRI data to standard BIDS format.
# Author: Krisanne Litinas

import os
import pydicom
import subprocess
import sys
import tarfile
import shutil
import json
import nipype.interfaces.fsl as fsl
import nibabel as nb
import datetime
import argparse
import glob

# Quick argument parser wrapper
def parse_args():
    parser = argparse.ArgumentParser(description='Convert MIDS data to standard BIDS format.')
    parser.add_argument('SUBJDIR', metavar='SubjectDir', type=str, help='Subject directory to be converted')
    parser.add_argument('--subj_id',help='Output subject ID to use for BIDS naming (default SubjectDir)')
    parser.add_argument('--build_cfg',action='store_true',help='Build dcm2bids config file')
    parser.add_argument('--convert',action='store_true',help='Run the converter')
    parser.add_argument('--d2b_cfg_file',help='Path of config file for dcm2bids to use. If none specified, generates one.')
    parser.add_argument('--bval_file',help='Path to .bval file to use for DTI dataset. Must have cooresponding .bvec')
    parser.add_argument('--exception-file',help='Path to file containing variables to override defaults')
    args = parser.parse_args()
    return(args)

# Set local environment
def set_env():
    os.environ['FSLOUTPUTTYPE'] = 'NIFTI_GZ'

    # Set this as wherever the BIDS data should land.
    global ROOTBIDSDIR
    ROOTBIDSDIR = '/export/klitinas/bids/root'

    if not os.path.isdir(ROOTBIDSDIR) or not os.path.isfile(os.path.join(ROOTBIDSDIR,'participants.tsv')):
        print("\nSetting up BIDS scaffolding in {}\n".format(ROOTBIDSDIR))
        make_bids_dirs(ROOTBIDSDIR)

        SCRIPTPATH = os.path.dirname(os.path.abspath(__file__))

        # Copy over config directory if needed.
        if not os.path.isdir(os.path.join(ROOTBIDSDIR,'code','config')):
            shutil.copytree(os.path.join(SCRIPTPATH,'config'),os.path.join(ROOTBIDSDIR,'code','config'))


        with open(os.path.join(ROOTBIDSDIR,'.bidsignore'),'w') as f:
            f.write('scratch/')
            #f.write('\nconfig/')

# Try to get run type
def runtype(FILEROOT):
    NII = FILEROOT + '.nii.gz'
    JSONFILE = FILEROOT + '.json'
    JSON_DAT = read_json(JSONFILE)
    DAT = nii_img(NII)
    HDR = nii_hdr(NII)

    RUNTYPE = 'anat'

    # DTI
    if os.path.isfile(FILEROOT + '.bval'):
        if DAT.shape[3] < 10:
            return 'dti_fm'
        else:
            return 'dti'
    else:
        # If not DTI but is hyperband, assume func
        if 'HYPERBAND' in JSON_DAT['ScanOptions']:
            return 'func_mb'

        if 'func_' in JSON_DAT['SeriesDescription']:
            return 'func_epi'

        if 'FM_' in JSON_DAT['SeriesDescription']:
            return 'func_fm'


    return RUNTYPE # default

# Get nii data
def nii_img(NIIFILE):
    return nb.load(NIIFILE)

# Nifti header
def nii_hdr(NIIFILE):
    opener = nb.openers.ImageOpener(NIIFILE)
    hdr = nb.Nifti1Header.from_fileobj(opener, check=False)
    return hdr

# Write out json file given dictionary
def write_json(FILENAME,JSONDICT):
    with open(FILENAME,'w') as f:
        json.dump(JSONDICT, f,indent=4)

# Returns json object from .json file
def read_json(FILE):
    with open(FILE,'r') as f:
        JSON = json.load(f)

    return JSON

# Gets directory which stores different aspects of configuration
def config_loc():
    return os.path.join(ROOTBIDSDIR,'code','config')


# sidecarChanges part of config for a func task
def sidecar_changes_func(HELPER_STR):

    HELPER_JSON_FILE = HELPER_STR + '.json'
    HELPER_NIIGZ = HELPER_STR + '.nii.gz'
    JSON_HELP = read_json(HELPER_JSON_FILE)

    RUNTYPE = runtype(HELPER_STR)

    hdr = nii_hdr(HELPER_NIIGZ)
    NUMSLICES = hdr.get_n_slices()

    if RUNTYPE == 'func_mb':

        sidecarChanges = {
                    'PulseSequenceType': "Multiband gradient echo EPI", \
                    'TaskName': JSON_HELP['SeriesDescription'], \
                    'SliceTiming': mb_slicetimes(JSON_HELP['MultibandAccelerationFactor'],NUMSLICES),
                    'PhaseEncodingDirection': "j-"
                        }

    else:
        sidecarChanges = {
                    'PulseSequenceType': "EPI", \
                    'TaskName': JSON_HELP['SeriesDescription']
                    #'TotalReadoutTime': TERT
        }

    return sidecarChanges

# Stub to get slice time list for mb
# Currently just has mb3/17sl, mb6/60sl
def mb_slicetimes(MBFACTOR,NUMSLICES):
    ST_LUT = os.path.join(config_loc(),'slice_timing.json')
    ST_JSON = read_json(ST_LUT)

    JSONFIELD = 'mb{:02d}_{}sl'.format(MBFACTOR,NUMSLICES)
    SLICETIMES = ST_JSON[JSONFIELD]
    return SLICETIMES

# Strips volumes from multiband .nii (default, strips first 12)
def mb_stripvols(NII,OUT_NII=None, STARTVOL=12):
    # Default: overwrite input file
    if not OUT_NII:
        OUT_NII = NII

    # TODO: figure out weird error message (if new file, do a touch first to create)
    fslroi = fsl.ExtractROI(in_file=NII, roi_file=OUT_NII, t_min=STARTVOL,t_size=-1)
    fslroi.run()

# Make some needed tweaks in fmap bidsdir
def fix_fmap(FMAPDIR = os.getcwd()):
    PDIR = os.getcwd()
    os.chdir(FMAPDIR)

    # TODO figure out DTI field map volumes
    for FILE in sorted(os.listdir('.')):
        BASE = os.path.splitext(FILE)[0]
        if FILE.endswith('.json') and not os.path.isfile(FILE + '.bval'):
            modify_func_fm_run(BASE)
        if FILE.endswith('.bval') or FILE.endswith('.bvec'):
            os.remove(FILE)

    os.chdir(PDIR)

# Stub to create the pepolar set of fm, given root file
def modify_func_fm_run(FILEROOT):
    # Input is like sub-sub_run-01_epi.nii.gz (already in fmap dir)
    NII = FILEROOT + '.nii.gz'
    JSONFILE = FILEROOT + '.json'
    JSON_DAT = read_json(JSONFILE)
    TASKNAME = JSON_DAT['SeriesDescription'].split('_')[1]

    # Splits into AP/PA sets
    # sub-<label>[_ses-<label>][_acq-<label>][_ce-<label>]_dir-<label>[_run-<index>]_epi.json
    LROOT = FILEROOT.split('_')
    APNAME = LROOT[0]+'_dir-AP_'+ '_'.join(LROOT[1:])
    PANAME = LROOT[0]+'_dir-PA_'+ '_'.join(LROOT[1:])


    # Volumes 1/3 are AP, 2/4 PA
    SPLITTER = fsl.Split(in_file = NII,dimension = 't',out_base_name = 'tmp')
    SPLITTER.run()

    MERGE_AP = fsl.Merge(in_files = ['tmp0000.nii.gz','tmp0002.nii.gz'],dimension='t',merged_file=APNAME+'.nii.gz')
    MERGE_AP.run()

    MERGE_PA = fsl.Merge(in_files = ['tmp0001.nii.gz','tmp0003.nii.gz'],dimension='t',merged_file=PANAME+'.nii.gz')
    MERGE_PA.run()

    os.remove(NII)
    os.remove(JSONFILE)
    for FILE in os.listdir('.'):
        if FILE.startswith('tmp0'):
            os.remove(FILE)

    # TODO: should we flip 1st volume and reorient?


    # Find which func data to use it for (right now assumes series description matches, may need to be more open).
    # Kludge: if stub of IntendedFor is there, populate with all runs
    os.chdir('..')

    INTENDEDFOR = []
    if 'IntendedFor' in JSON_DAT:
        STUB_TASKS = JSON_DAT['IntendedFor']
        for STUB in STUB_TASKS:
            INTEND_TASK = glob.glob('func/*task-{}_*bold.nii.gz'.format(STUB))
            INTENDEDFOR += INTEND_TASK
    else:
        INTENDEDFOR = glob.glob('func/*task-{}_*bold.nii.gz'.format(TASKNAME))

    JSON_DAT['IntendedFor'] = sorted(INTENDEDFOR)
    os.chdir('fmap')

    write_json(APNAME + '.json',JSON_DAT)
    write_json(PANAME + '.json',JSON_DAT)

# Quick wrapper to dcm2bids_scaffold util
def make_bids_dirs(OUTDIR=os.getcwd()):

    TMP = subprocess.Popen(['dcm2bids_scaffold', '-o', OUTDIR], stdout = subprocess.PIPE)
    OUTPUT = str(TMP.communicate())

    # README can't be empty
    with open(os.path.join(OUTDIR,'README'),'w') as f:
        f.write('Converted with mids2bids')

    # dataset_description.json funding needs to be list, not str
    DESCFILE = os.path.join(OUTDIR,'dataset_description.json')
    DESC = read_json(DESCFILE)
    DESC['Funding'] = [""]
    write_json(DESCFILE,DESC)

# Stub to find the appropriate config file to use
def dcm2bids_study_config(SUBJDIR=None):
    # Get based on protocol name
    if not SUBJDIR:
        SUBJDIR = os.getcwd()

    CONFIGDIR = os.path.join(config_loc(),'.dcm2bids_config')
    PAT = '{}/sourcedata/dicom/*/i*.1'.format(SUBJDIR)
    DCM = glob.glob(PAT)[0]
    d = pydicom.read_file(DCM)
    CFGNAME = os.path.join(config_loc(),'.dcm2bids_config',d.ProtocolName.replace(' ','_').replace('/','.')+'.json')
    return CFGNAME

# Wrapper for dcm2bids call
def run_dcm2bids(SUBJECT, CONFIGFILE):

    TMP = subprocess.Popen(['dcm2bids', '-d', 'sourcedata/dicom/', '-p', SUBJECT, '-c', CONFIGFILE], stdout = subprocess.PIPE)
    OUTPUT = str(TMP.communicate())
    return OUTPUT

# Massage dcm2bids outputs so they pass validation
def finetune(PARENTDIR,SUBJID):

    # Deal with field map data
    if os.path.isdir(os.path.join(PARENTDIR,'sub-' + SUBJID,'fmap')):
        fix_fmap(os.path.join(PARENTDIR,'sub-' + SUBJID,'fmap'))

# Copy over dicoms and (eventually) pfiles and extract
def grab_source_files(ORIGDIR,OUTDIR=os.getcwd()):

    # DICOMS
    DCMTGZ = os.path.join(ORIGDIR,'dicom.tgz')
    if os.path.isfile(DCMTGZ):
        with tarfile.open(DCMTGZ,"r:gz") as tar:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar, path=OUTDIR)
    # Dicoms, if already extracted.
    elif os.path.isdir(os.path.join(ORIGDIR,'dicom')):
        shutil.copytree(os.path.join(ORIGDIR,'dicom'),OUTDIR)

    else:
        print("No dicoms found in {}!".format(ORIGDIR))

    # TODO: PFILES
    #if os.path.isdir(os.path.join(ORIGDIR,'raw','pfiles')):
    #    shutil.copytree(os.path.join(ORIGDIR,'raw','pfiles'),OUTDIR)

def run_dcm2bids_helper():
    TMP = subprocess.Popen(['dcm2bids_helper', '-d', 'sourcedata/dicom'], stdout = subprocess.PIPE)
    OUTPUT = str(TMP.communicate())

# Stub to build as-best-as-possible config based on dcm2bids_helper etc
def build_config(SUBJDIR=None,EXTRAFILE=None):

    if not SUBJDIR:
        SUBJDIR = os.getcwd()

    if EXTRAFILE:
        EXTRADAT = read_json(EXTRAFILE)
    else:
        EXTRADAT = {}

    # make directory structure
    BIDSDIR = os.path.join(ROOTBIDSDIR,'scratch',os.path.basename(SUBJDIR))

    print("Making directory structure in {}".format(BIDSDIR))

    if not os.path.isdir(BIDSDIR):
        os.makedirs(BIDSDIR)

    if not os.path.isdir(os.path.join(BIDSDIR,'sourcedata')):
        make_bids_dirs(BIDSDIR)

    # mv dicoms into sourcedata
    if not os.path.isdir(os.path.join(BIDSDIR,'sourcedata','dicom')):
        print("Copying raw data from {} to {}".format(SUBJDIR,os.path.join(BIDSDIR,'sourcedata')))
        grab_source_files(SUBJDIR,os.path.join(BIDSDIR,'sourcedata'))

    # 1. run dcm2bids_helper
    print("Changing directories to {}".format(BIDSDIR))
    os.chdir(BIDSDIR)
    if not os.path.isdir(os.path.join(BIDSDIR,'tmp_dcm2bids/helper')):
        print("Running dcm2bids_helper script.")
        run_dcm2bids_helper()

    print("Generating config file")

    # 2. loop through the resulting json data to write the config

    # Construct output file name
    CFGNAME = dcm2bids_study_config()

    # First, stick generic stubs for anat, fieldmaps, and dti
    CONFIG_OUT = read_json(os.path.join(config_loc(),'anat.json'))
    FMCFG = read_json(os.path.join(config_loc(),'fieldmap.json'))
    DTICFG = read_json(os.path.join(config_loc(),'dti.json'))

    CONFIG_OUT['descriptions'] += FMCFG['descriptions']
    CONFIG_OUT['descriptions'] += DTICFG['descriptions']

    # Loop through dcmhelper-generated json files
    os.chdir('tmp_dcm2bids/helper')
    for JSONFILE in glob.glob('*json'):
        BASEFILE = os.path.splitext(JSONFILE)[0]
        RUNTYPE = runtype(BASEFILE)

        JSDAT = read_json(JSONFILE)
        SEDESC = JSDAT['SeriesDescription']

        # func, along with any specified exceptions
        if RUNTYPE.startswith('func_') and RUNTYPE != "func_fm":

            RUNCFG = func_task_config(BASEFILE)
            if SEDESC in EXTRADAT.keys():
                if 'sidecarChanges' not in RUNCFG.keys():
                    RUNCFG['sidecarChanges'] = {}

                for ITEM in EXTRADAT[SEDESC]:
                    RUNCFG['sidecarChanges'][ITEM] = EXTRADAT[SEDESC][ITEM]


            if RUNCFG not in CONFIG_OUT['descriptions']:
                print(SEDESC)
                CONFIG_OUT['descriptions'].append(RUNCFG)

    # Write out the file, but don't overwrite if exists
    if os.path.isfile(CFGNAME):
        print("File {} already exists, moving to backup location.".format(CFGNAME))
        DATESTR = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        BASENAME = os.path.splitext(os.path.basename(CFGNAME))[0]
        BACKUPFILE = os.path.join(os.path.dirname(CFGNAME),'previous',DATESTR+'_'+BASENAME+'.json')
        shutil.move(CFGNAME,BACKUPFILE)

    write_json(CFGNAME,CONFIG_OUT)

    print("Config generated: {}".format(CFGNAME))

# builds config stub for one func task (mb)
def func_task_config(HELPER_STR):
    JSON = read_json(HELPER_STR + '.json')
    TASK = JSON['SeriesDescription'].split('_')[1]
    CONFIGOBJ = {
            "dataType": "func",
            "modalityLabel": "bold",
            "customLabels": "task-{}".format(TASK),
            "criteria": {
                "SeriesDescription": JSON['SeriesDescription']
              },
              "sidecarChanges": sidecar_changes_func(HELPER_STR)
    }

    return CONFIGOBJ

# Pad with extra zeros for B0 volumes (needed for multiband)
def pad_bvec_bval_file(BVFILE):
    NIIGZ = os.path.splitext(BVFILE)[0] + '.nii.gz'
    NIIHDR = nii_hdr(NIIGZ)
    NVOLS = NIIHDR.get_data_shape()[3]

    # Pad each line with needed number of zeros
    with open(BVFILE,'r+') as f:
        DAT = f.readlines()

        f.seek(0)
        for LINE in DAT:
            VALS = LINE.strip()
            NVALSORIG = len(VALS.split(' '))
            if NVALSORIG < NVOLS:
                DIFF = NVOLS - NVALSORIG
                PREFIX = '0 ' * DIFF
                OUTVAL = PREFIX + VALS + "\n"
                f.write(OUTVAL)


if __name__ == "__main__":
    myargs = parse_args()

    # Set a few variables
    set_env()

    # Build config file if directed
    if myargs.build_cfg:
        build_config(myargs.SUBJDIR,myargs.exception_file)

    # If output subj name not specified, make it the subject dir name (minus any '_')
    if myargs.subj_id:
        OUTSUBJNAME = myargs.subj_id
    else:
        OUTSUBJNAME = os.path.basename(myargs.SUBJDIR).replace('_','')

    # Run conversion if directed
    # TODO (maybe): generate bids dirs (SUBJ_BIDS) + sourcedata subdirs if not already there
    if myargs.convert:

        # TODO: repeats in build_cfg, should maybe make own fcn
        WORKDIR = os.path.join(ROOTBIDSDIR,'scratch',os.path.basename(myargs.SUBJDIR))
        if not os.path.isdir(WORKDIR):
            os.makedirs(WORKDIR)

            print("Making directory structure in {}".format(WORKDIR))

            if not os.path.isdir(os.path.join(WORKDIR,'sourcedata')):
                make_bids_dirs(WORKDIR)

            # mv dicoms into sourcedata
            if not os.path.isdir(os.path.join(WORKDIR,'sourcedata','dicom')):
                print("Copying raw data from {} to {}".format(myargs.SUBJDIR,os.path.join(WORKDIR,'sourcedata')))
                grab_source_files(myargs.SUBJDIR,os.path.join(WORKDIR,'sourcedata'))

        os.chdir(WORKDIR)

        # If config file specified, use that. otherwise look for it and if not found, prompt to make one.
        if myargs.d2b_cfg_file:
            CFGFILE =  myargs.d2b_cfg_file
        else:
            CFGFILE = dcm2bids_study_config()
            if not os.path.isfile(CFGFILE):
                BUILD_CFG = input("\nExisting config file {} not found - build it? (Y/N) \n".format(CFGFILE))
                if BUILD_CFG.lower() == 'y':
                    build_config(myargs.SUBJDIR,myargs.exception_file)
                    print("Check it over before running converter.")
                    sys.exit()

                elif BUILD_CFG.lower() == 'n':
                    print("\nOkay, exiting.")
                    sys.exit()

        # Run dcm2bids conversion
        run_dcm2bids(OUTSUBJNAME, CFGFILE)

        # maybe put this after dti stuff
        finetune(WORKDIR,OUTSUBJNAME)

        # Replace bval/bvec auto-generated by dcm2niix
        if myargs.bval_file:
            BVECFILE = myargs.bval_file.replace('.bval','.bvec')

            DTIDIR = os.path.join(WORKDIR,'sub-'+OUTSUBJNAME,'dwi')
            if os.path.isdir(DTIDIR):
                BVALNAME = glob.glob('{}/*bval'.format(DTIDIR))
                BVECNAME = glob.glob('{}/*bvec'.format(DTIDIR))

                for MYBVAL in BVALNAME:
                    MYBVEC = MYBVAL.replace('.bval','.bvec')
                    os.remove(MYBVAL)
                    shutil.copyfile(myargs.bval_file,MYBVAL)
                    shutil.copyfile(BVECFILE,MYBVEC)

                    pad_bvec_bval_file(MYBVAL)
                    pad_bvec_bval_file(MYBVEC)

        # Move results to main subject directory
        shutil.move(os.path.join(WORKDIR,'sub-'+OUTSUBJNAME),ROOTBIDSDIR)

        # Copy source data and processed results into the BIDS dir structure
