# mids2bids
Work in progress to convert UMICH MRI data to BIDS format.

Attempts to auto-generate a dcm2bids configuration file which can then be used to do the actual file conversion.

## Setup
### Software requirements
* dcm2niix (tested with v.1.0.20200331)
* dcm2bids (https://cbedetti.github.io/Dcm2Bids/, tested with v. 2.1.4)
* python3 with these packages:
    - nipype
    - nibabel
    - pydicom

### Set local environment
*  Download `mids2bids.py` and the `config` directory to the same parent folder.  `config` will eventually be copied into wherever the BIDS subject directory is set.
*  Modify line in set_env function (line #38 in current version) of `mids2bids.py` point to directory where the BIDSified datasets will land.

## Usage
Display help text
```
mids2bids.py -h  
```
General Syntax
```
mids2bids.py <SUBJDIR_TO_BE_CONVERTED> <opts>
```

### Making a configuration file
A main component of the conversion is based on dcm2bids, which requires a configuration file to sort the data. This config file is generally protocol-specific.  Assuming no changes across subjects, you should ideally only need to do this once, and then use that config file for converting all data. 

To build a config file, point to a representative subject for the study and use the `--build_cfg` flag:
```
mids2bids.py <SUBJDIR_TO_BE_CONVERTED> --build_cfg
```
A .json file will be created in the [bidsroot]/code/config/.dcm2bids_config directory, named by protocol name.  If the file already exists, will backup the old one.

The config file is made by using the `dcm2bids_helper` utility along with some other customizations specific to UMICH data.  Additionally the `exception-file` argument may also be used, specifying a json file location that contains known overrides to what the config would usually contain that cannot be generated automatically.  For instance, `TotalReadoutTime` values of dcm2niix-generated .json sidecars are unreliable for GE data.  They also cannot be pulled from DICOM, so must be explicitly listed.

For an example of the `exception-file` layout, see config/exceptions.json.

### Running the conversion
To run the actual conversion, use the `--convert` flag:
```
mids2bids.py <SUBJDIR_TO_BE_CONVERTED> --convert
```
If the `--d2b_cfg_file` optional argument is not specified, looks for a config file to use in [bidsroot]/code/config/.dcm2bids_config/.  If not found for this protocol, will prompt user to generate one.

### Other optional arguments
* **subj-id**: Subject ID to use for BIDS outputs.  Resulting BIDS directory will be called 'sub-[subj-id]'. Default is the basename of the input subject directory.
* **bval_file**:Path of bval file to use instead of dcm2niix-generated one, which is incorrect for multiband data. Requires a corresponding .bvec file.
* **d2b_cfg_file**: Path of dcm2bids config file to use for conversion. If not specified, will try to automatically find it, and if not detected will prompt to build it.
* **exception-file**: used for building the dcm2bids config file. json file containing key/values that should be used instead of any defaults. Key names are series descriptions.  An example is located at config/exceptions.json


## Future work / current questions / things to note
* Add pfile compatibility, currently just DICOM-friendly.
* Need to fully verify PEPOLAR, timing parameters for functional and DTI data. Current work is based mainly on ABCD data.
* DTI field map data is probably/certainly incorrectly split between AP/PA PhaseEncodingDirections
* Functional field maps are mapped to the corresponding fMRI data (the "IntendedFor" field in the .json). This may not always match exactly.
* Final output BIDS directory do not include our processed outputs or source DICOMS/pfiles. Could easily do a copy from original (processed data into BIDS derivitives/ folder, DICOMS/pfiles into sourcedata/ would work), but unsure how best we want to handle that.  Maybe have that as an option.
* Might sometimes run into problems if you try to repeat the conversion for a given subject without removing outputs from the first attempt beforehand.
* dcm2niix, when called from dcm2bids, sometimes hangs and takes some time to finish.  I believe this mainly happens with high disk I/O and/or low available disk space.
