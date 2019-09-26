#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf                                             Path to vcf file

Optional parameters:

GENERAL

  VCF_2_FASTA
    -h|--help                                            Shows help
    -f|--flank                                           Flanking size to design primer [$FLANK]
    -o|--offset                                          Offset, one-sided amplicon size [$OFFSET]
    -m|--mask                                            Mask the sequence [$MASK]

    -vfs|--vcf_fasta_script                              Path to vcf_to_fasta.py script [$VCF_FASTA_SCRIPT]
    -venv|--venv                                         Path to virtual environment [$VENV]

  PRIMER_DESIGN
    -pt|--pcr_type                                       PCR type [$PCR_TYPE]
    -tp|--tilling_params                                 Tilling parameters [$TILLING_PARAMS]
    -psr|--psr                                           PSR [$PSR]
    -mp|-mispriming                                      Mispriming file [$MISPRIMING]

    -bd|--bindir                                         Bindir [$BINDIR]
    -gp|--guix_profile                                   GUIX profile [$GUIX_PROFILE]
    -pc|--primer3_core                                   Primer3 core [$PRIMER3_CORE]

"
exit
}

POSITIONAL=()

# DEFAULTS
FLANK=200
OFFSET=0
MASK=false
VCF_FASTA_SCRIPT='/hpc/pmc_vanboxtel/tools/PrimerDesign/vcf_to_fasta.py'
VENV='/hpc/pmc_vanboxtel/tools/PrimerDesign/venv_3.6/bin/activate'

PCR_TYPE='single'
TILLING_PARAMS=' '
PSR='30-($FLANK+30)'
BINDIR='/hpc/pmc_vanboxtel/tools/PrimerDesign/primer3/primers'
GUIX_PROFILE='/hpc/pmc_vanboxtel/tools/PrimerDesign/primer3/emboss/.guix-profile'
PRIMER3_CORE='/hpc/pmc_vanboxtel/tools/PrimerDesign/primer3/primer3/src/primer3_core'
MISPRIMING='/hpc/pmc_vanboxtel/tools/PrimerDesign/primer3/repbase/current/empty.ref'
#MISPRIMING=' /hpc/pmc_vanboxtel/tools/PrimerDesign/primer3/repbase/current/humrep.ref'

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
# REQUIRED OPTIONS
    -v|--vcf)
    VCF="$2"
    shift # past argument
    shift # past value
    ;;
# GENERAL OPTIONS
    -h|--help)
    usage
    shift # past argument
    ;;
# VCF_2_FASTA OPTIONS
    -f|--flank)
    FLANK="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--offset)
    OFFSET="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--mask)
    MASK=true
    shift # past argument
    ;;
    -vfs|--vcf_fasta_script)
    VCF_FASTA_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -venv|--venv)
    VENV="$2"
    shift # past argument
    shift # past value
    ;;
# PRIMER DESIGN OPTIONS
    -pt|--pcr_type)
    PCR_TYPE="$2"
    shift # past argument
    shift # past value
    ;;
    -tp|--tilling_params)
    TILLING_PARAMS="$2"
    shift # past argument
    shift # past value
    ;;
    -psr|--psr)
    PSR="$2"
    shift # past argument
    shift # past value
    ;;
    -mp|--mispriming)
    MISPRIMING="$2"
    shift # past argument
    shift # past value
    ;;
    -bd|--bindir)
    BINDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -gp|--guix_profile)
    GUIX_PROFILE="$2"
    shift # past argument
    shift # past value
    ;;
    -pc|--primer3_core)
    PRIMER3_CORE="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
if [ -z $VCF ]; then
  echo "Missing -v|--vcf parameter"
  usage
fi

if [ ! -f $VCF ]; then
  echo "VCF file \"$VCF\" does not exists."
  exit
fi

prepare() {
  if [[ $PSR =~ "FLANK" ]]; then
    PSR='30-'$((FLANK+30))
  fi
  if [ -d $PWD/primer3 ]; then
    rm -rf $PWD/primer3
  else
    mkdir $PWD/primer3
  fi
  FASTA=${VCF/.vcf/.fasta}
  FASTA=$(basename $FASTA)
  FASTA=$PWD/$FASTA
  LOG=${VCF/.vcf/_primer3.log}
  ERR=${VCF/.vcf/_primer3.err}

}

vcf_2_fasta() {
  . $VENV

  if [ $MASK = true ]; then
    python $VCF_FASTA_SCRIPT -o $OFFSET -f $FLANK -m $VCF > $FASTA
  else
    python $VCF_FASTA_SCRIPT -o $OFFSET -f $FLANK $VCF > $FASTA
  fi
}

primer_design() {
cd $PWD/primer3

guixr load-profile $GUIX_PROFILE --<<EOF

export EMBOSS_PRIMER3_CORE=$PRIMER3_CORE
$BINDIR/primerBATCH1 $MISPRIMING $PCR_TYPE $PSR $TILLING_PARAMS <$FASTA 1>$LOG 2>$ERR
EOF

cd $PWD
}

prepare
vcf_2_fasta
primer_design

exit
