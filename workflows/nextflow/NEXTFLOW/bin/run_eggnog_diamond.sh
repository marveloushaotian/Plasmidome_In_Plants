#!/bin/bash
# Automatic parallel EggNOG diamond search script

set -e

# Parse arguments
MAX_PARALLEL=10  # Number of parallel tasks, adjust based on available resources
CHUNK_DIR=""
DATABASE_DIR=""
CPU_PER_TASK=20

while [[ $# -gt 0 ]]; do
    case $1 in
        --chunk_dir)
            CHUNK_DIR="$2"
            shift 2
            ;;
        --database)
            DATABASE_DIR="$2"
            shift 2
            ;;
        --max_parallel)
            MAX_PARALLEL="$2"
            shift 2
            ;;
        --cpu)
            CPU_PER_TASK="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ -z "$CHUNK_DIR" ] || [ -z "$DATABASE_DIR" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --chunk_dir <dir> --database <dir> [--max_parallel <n>] [--cpu <n>]"
    exit 1
fi

echo "Starting EggNOG diamond search..."
echo "Chunk directory: $CHUNK_DIR"
echo "Database directory: $DATABASE_DIR"
echo "Max parallel tasks: $MAX_PARALLEL"
echo "CPU per task: $CPU_PER_TASK"

cd "$CHUNK_DIR"

# Process each chunk file
for chunk_file in *.part_*.faa; do
    # Check if file exists (in case glob matches nothing)
    [ -e "$chunk_file" ] || continue
    
    # Control concurrent jobs
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
        sleep 5  # Wait 5 seconds before checking again
    done
    
    echo "Processing: $chunk_file"
    emapper.py -m diamond --no_annot --no_file_comments --cpu $CPU_PER_TASK \
        -i "$chunk_file" -o "$chunk_file" \
        --data_dir "$DATABASE_DIR" -d none --itype proteins --override \
        1>"${chunk_file}.log" 2>&1 &
done

# Wait for all tasks to complete
wait
echo "All diamond searches completed!"

