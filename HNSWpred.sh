# usage function
usage() {
    echo "Usage: $0 <data_dir> <top_k> <user_file> <user_output_file>"
    exit 1
}

if [[ -d $1 ]]; then
    re='^[0-9]+$'
    if [[ $2 =~ $re ]] ; then
        if [[ -f $3 ]]; then
            make run data_dir=$1 top_k=$2 user=$3 out_file=$4
        else
            usage
        fi
    else
        usage
    fi
else
    usage
fi