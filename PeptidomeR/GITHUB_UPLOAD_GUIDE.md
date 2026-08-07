# Upload PeptidomeR to GitHub

Repository: [indrani2611/PeptidomeR](https://github.com/indrani2611/PeptidomeR)

The safest approach is to clone the existing GitHub repository and copy the prepared project into that clone. This keeps any commits already present on GitHub and avoids overwriting the remote history.

## 1. Install and configure Git

Open Terminal on the Mac and check Git:

```bash
git --version
```

If Git is not available, macOS should offer to install the Command Line Tools. Configure the name and email used in commits if this has not been done previously:

```bash
git config --global user.name "Indrani Bera"
git config --global user.email "YOUR_GITHUB_EMAIL"
```

Replace `YOUR_GITHUB_EMAIL` with the email used for GitHub, or with the GitHub-provided private no-reply email.

## 2. Download and unpack the prepared project

Download `PeptidomeR_GitHub_ready.zip` from this conversation and double-click it in Finder. This creates a `PeptidomeR_GitHub_ready` folder, normally in `Downloads`.

## 3. Clone the GitHub repository

In Terminal:

```bash
cd ~/Desktop
git clone https://github.com/indrani2611/PeptidomeR.git
cd PeptidomeR
```

If the repository is empty, Git may display a warning. That is normal.

## 4. Copy the prepared files into the clone

```bash
rsync -av --exclude='.git' ~/Downloads/PeptidomeR_GitHub_ready/ ./
```

If the extracted folder is in a different location, drag the folder from Finder into Terminal after typing `rsync -av --exclude='.git' `, then add a space and `./`.

## 5. Review what Git will upload

```bash
git status
```

You should see `app.R`, `R/`, `www/`, `docs/`, `Dockerfile`, `renv.lock`, `README.md`, and the configuration files. `.DS_Store`, `__MACOSX`, local R state, and ZIP archives should not be listed.

## 6. Commit the application

```bash
git add .
git status
git commit -m "Add PeptidomeR Shiny application and documentation"
```

Read the second `git status` output before committing to ensure no private datasets, passwords, tokens, or `.Renviron` files are staged.

## 7. Push to GitHub

```bash
git branch -M main
git push -u origin main
```

Complete GitHub authentication in the browser if prompted. Do not paste a password or access token into source files.

## 8. Verify the repository

Open [https://github.com/indrani2611/PeptidomeR](https://github.com/indrani2611/PeptidomeR) and confirm that:

- the logo appears in the main README;
- all six files are present under `R/`;
- the module guides open from the README table;
- `renv.lock` and `Dockerfile` are present; and
- the GitHub Actions workflow begins a Docker build under the **Actions** tab.

## Later updates

For subsequent changes, work inside the cloned `~/Desktop/PeptidomeR` folder:

```bash
cd ~/Desktop/PeptidomeR
git pull --rebase origin main
git status
git add .
git commit -m "Describe the update"
git push origin main
```

Replace `Describe the update` with a short summary such as `Fix multivariate heatmap display`.

## Common problems

### Push rejected because the remote contains newer work

```bash
git pull --rebase origin main
git push origin main
```

If Git reports a conflict, stop and inspect the named files before continuing. Do not use a force push unless you intentionally want to replace remote history.

### A file larger than 100 MB is rejected

Do not commit analysis datasets, generated ZIP results, local package libraries, or container images. The source repository should contain code, documentation, the logo, and dependency files.

### Docker test says port 3838 is already allocated

Use another local port:

```bash
docker run --rm -p 3839:3838 peptidomer:latest
```

Open [http://localhost:3839](http://localhost:3839).

## Before announcing the repository

Add the final paper citation and an institutionally approved open-source licence. The licence should be chosen with all authors and UCD requirements in mind.

