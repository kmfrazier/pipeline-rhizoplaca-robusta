#!/usr/bin/env bash
set -e

# Configure Git identity
git config --global user.name "rachel-hunter-smhs"
git config --global user.email "rachel.hunter2001@gmail.com"

# Define paths
REPO_DIR="$HOME/pipeline-rhizoplaca-robusta"
WORKFLOW_SRC="/grphome/grp_lichenscapstone/blastx3-5/blastx_nitrogenfixation"
DEST_DIR="$REPO_DIR/nitrogen_fixation"
DB_SRC="/grphome/grp_lichenscapstone/combined_db_uniprot"
DB_PREFIX="combined_db"
CHUNK_SCRIPT_SRC="/grphome/grp_lichenscapstone/chunkSmallBatch.slurm"

# 1. Clone your fork if not already present
if [ ! -d "$REPO_DIR" ]; then
  cd "$HOME"
  git clone git@github.com:rachel-hunter-smhs/pipeline-rhizoplaca-robusta.git
  # Ensure remote uses SSH URL
git -C "$HOME/pipeline-rhizoplaca-robusta" remote set-url origin git@github.com:rachel-hunter-smhs/pipeline-rhizoplaca-robusta.git
fi

# 2. Create destination directory inside your repo
mkdir -p "$DEST_DIR"

# 3. Copy specific workflow scripts into the destination
if [ -f "$WORKFLOW_SRC/setup.sh" ]; then
  cp "$WORKFLOW_SRC/setup.sh" "$DEST_DIR/"
elif [ -f "/grphome/grp_lichenscapstone/setup.sh" ]; then
  cp "/grphome/grp_lichenscapstone/setup.sh" "$DEST_DIR/"
fi

cp "$WORKFLOW_SRC/organismInfo.py" "$DEST_DIR/"
cp "$WORKFLOW_SRC/process_blastx.py" "$DEST_DIR/" 2>/dev/null || true
cp "$WORKFLOW_SRC/graph.py" "$DEST_DIR/"

if [ -f "$CHUNK_SCRIPT_SRC" ]; then
  cp "$CHUNK_SCRIPT_SRC" "$DEST_DIR/"
fi

# 4. Recursively copy all other files and subdirectories, excluding logs
rsync -av --exclude='*.txt' --exclude='*.out' --exclude='*.err' "$WORKFLOW_SRC"/ "$DEST_DIR"/

# 5. Copy BLAST database files into a subdirectory
mkdir -p "$DEST_DIR/database"
cp "$DB_SRC/${DB_PREFIX}"* "$DEST_DIR/database/"

# 6. Commit and push to your fork
cd "$REPO_DIR"
git add nitrogen_fixation
# Only commit if there are staged changes
git diff --cached --quiet || git commit -m "Add nitrogen fixation workflow and database files"

git push origin main || git push origin HEAD:nitrogen_fixation_branch

echo "Deployment complete: files in $DEST_DIR"

