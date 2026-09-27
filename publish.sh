#!/usr/bin/env bash
# Publish this folder to GitHub and turn on GitHub Pages (website served from docs/).
# Run from inside the omics-ml-explorer folder:   bash publish.sh
set -e
REPO="omics-ml-explorer"
OWNER=$(gh api user -q .login)

if [ ! -d .git ]; then
  git init -q -b main
fi
git add -A
git commit -qm "Omics ML Explorer: Shiny app, dummy data and project website" || true

if gh repo view "$OWNER/$REPO" >/dev/null 2>&1; then
  git remote get-url origin >/dev/null 2>&1 || git remote add origin "https://github.com/$OWNER/$REPO.git"
  git push -u origin main
else
  gh repo create "$REPO" --public --source=. --push \
    --description "R Shiny app: protease classifiers (prostate cancer), ML feature selection (PDAC), longitudinal biomarkers (uveal melanoma)"
fi

# GitHub Pages from the docs/ folder of main
gh api -X POST "repos/$OWNER/$REPO/pages" \
  -f "source[branch]=main" -f "source[path]=/docs" >/dev/null 2>&1 \
|| gh api -X PUT "repos/$OWNER/$REPO/pages" \
  -f "source[branch]=main" -f "source[path]=/docs" >/dev/null 2>&1 || true

echo
echo "Repository: https://github.com/$OWNER/$REPO"
echo "Website:    https://$OWNER.github.io/$REPO/   (live in 1-3 minutes)"
echo "Run app:    in R ->  shiny::runGitHub(\"$REPO\", \"$OWNER\")"
