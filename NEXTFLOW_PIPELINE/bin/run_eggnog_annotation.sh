#!/bin/bash
# Parallel EggNOG annotation script with concurrency control

set -e

# =============Configuration Parameters=============
INPUT_DIR=""
OUTPUT_DIR=""
DATABASE_DIR=""
CPU_PER_TASK=10
TOTAL_CPU=120
MAX_CONCURRENT=$((TOTAL_CPU / CPU_PER_TASK))

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --database)
            DATABASE_DIR="$2"
            shift 2
            ;;
        --cpu_per_task)
            CPU_PER_TASK="$2"
            MAX_CONCURRENT=$((TOTAL_CPU / CPU_PER_TASK))
            shift 2
            ;;
        --total_cpu)
            TOTAL_CPU="$2"
            MAX_CONCURRENT=$((TOTAL_CPU / CPU_PER_TASK))
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$DATABASE_DIR" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --input_dir <dir> --output_dir <dir> --database <dir> [--cpu_per_task <n>] [--total_cpu <n>]"
    exit 1
fi

echo "=== Batch EggNOG-mapper Annotation Task ==="
echo "Total CPU: $TOTAL_CPU"
echo "CPU per task: $CPU_PER_TASK" 
echo "Max concurrent: $MAX_CONCURRENT"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Database directory: $DATABASE_DIR"
echo "=================================="

# Check if input folder exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Folder $INPUT_DIR does not exist!"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create file list
FILES=()
for file in "$INPUT_DIR"/*.emapper.seed_orthologs; do
    if [ -f "$file" ]; then
        FILES+=("$file")
    fi
done

if [ ${#FILES[@]} -eq 0 ]; then
    echo "Error: No .emapper.seed_orthologs files found"
    exit 1
fi

echo "Found ${#FILES[@]} files to process"
echo ""

# Store running task information
declare -A running_tasks  # Associative array: PID -> filename
current_file_index=0
completed_count=0
total_files=${#FILES[@]}

# Function: start single task
start_task() {
    local file="$1"
    local filename=$(basename "$file")
    local basename_no_ext=$(basename "$file" .emapper.seed_orthologs)
    local output_name="${basename_no_ext}_eggnog"
    local log_file="${OUTPUT_DIR}/annotation_${basename_no_ext}.log"
    
    echo "Starting task: $filename"
    echo "   Output prefix: $output_name"
    echo "   Log file: $log_file"
    
    # Start task
    emapper.py --annotate_hits_table "$file" \
               --no_file_comments \
               -o "${OUTPUT_DIR}/$output_name" \
               --cpu $CPU_PER_TASK \
               --data_dir "$DATABASE_DIR" \
               1>"$log_file" 2>&1 &
    
    local pid=$!
    running_tasks[$pid]="$filename"
    
    echo "   Started PID: $pid"
    echo "   Running tasks: ${#running_tasks[@]}/$MAX_CONCURRENT"
    echo ""
}

# Function: check and clean completed tasks
check_completed_tasks() {
    local completed_pids=()
    
    # Check each running task
    for pid in "${!running_tasks[@]}"; do
        if ! kill -0 "$pid" 2>/dev/null; then
            # Task completed
            completed_pids+=($pid)
            completed_count=$((completed_count + 1))
            echo "Task completed: ${running_tasks[$pid]} (PID: $pid)"
            echo "   Progress: $completed_count/$total_files"
        fi
    done
    
    # Remove completed tasks from running_tasks
    for pid in "${completed_pids[@]}"; do
        unset running_tasks[$pid]
    done
    
    echo "   Currently running: ${#running_tasks[@]}/$MAX_CONCURRENT"
    echo ""
}

# Function: wait until slot available
wait_for_available_slot() {
    while [ ${#running_tasks[@]} -ge $MAX_CONCURRENT ]; do
        echo "Max concurrency reached ($MAX_CONCURRENT), waiting for tasks to complete..."
        sleep 10  # Check every 10 seconds
        check_completed_tasks
    done
}

# Main loop: process all files
echo "Starting batch processing..."
echo ""

while [ $current_file_index -lt $total_files ]; do
    # Check and clean completed tasks
    check_completed_tasks
    
    # Wait for available slot
    wait_for_available_slot
    
    # Start new task
    file="${FILES[$current_file_index]}"
    start_task "$file"
    
    current_file_index=$((current_file_index + 1))
    
    # Show overall progress
    echo "Overall progress: Started $current_file_index/$total_files, Completed $completed_count/$total_files"
    echo "----------------------------------------"
done

# Wait for all remaining tasks to complete
echo "All tasks started, waiting for remaining tasks to complete..."
echo ""

while [ ${#running_tasks[@]} -gt 0 ]; do
    echo "Waiting for remaining ${#running_tasks[@]} tasks..."
    echo "   Running tasks:"
    for pid in "${!running_tasks[@]}"; do
        echo "     - PID $pid: ${running_tasks[$pid]}"
    done
    echo ""
    
    sleep 15  # Check every 15 seconds
    check_completed_tasks
done

# Final results
echo ""
echo "=============== Task Completed ==============="
echo "Total files processed: $total_files"
echo "Successfully completed: $completed_count"
echo ""

# Check generated files
echo "Checking result files:"
result_count=$(ls "$OUTPUT_DIR"/*_eggnog.emapper.annotations 2>/dev/null | wc -l)
if [ $result_count -gt 0 ]; then
    echo "Found $result_count annotation result files:"
    ls -1 "$OUTPUT_DIR"/*_eggnog.emapper.annotations 2>/dev/null | head -10
    if [ $result_count -gt 10 ]; then
        echo "   ... and $((result_count - 10)) more files"
    fi
else
    echo "No annotation result files found, please check logs"
fi

echo ""
echo "View detailed results:"
echo "   - All result files: ls -la $OUTPUT_DIR/*_eggnog.*"
echo "   - View log errors: grep -i error $OUTPUT_DIR/annotation_*.log"
echo "   - View completion status: grep -i 'done\|finished\|completed' $OUTPUT_DIR/annotation_*.log"
echo "==============================================="

