# Publish this repository on your GitHub account

The supplied ZIP is **ready for upload**, but it has not been published to
GitHub. Extract it first; do not upload just the ZIP if you want visitors to
browse the Python scripts and README directly.

## Option 1: GitHub Desktop (friendly for Windows)

1. Extract `NGS-Sequencing-QC-Analysis.zip` to a folder on your computer.
2. Open GitHub Desktop and sign in to your GitHub account.
3. Select **File > Add local repository** and browse to the extracted
   `NGS-Sequencing-QC-Analysis` folder. If GitHub Desktop says it is not a
   repository, select **create a repository here**.
4. Confirm the repository name `NGS-Sequencing-QC-Analysis`, commit the files,
   then select **Publish repository**.
5. Choose **Public** only if you want hiring managers to view it without a
   sign-in. Review the files before publishing; all included reads are synthetic.
6. The public URL will be in the form:
   `https://github.com/YOUR_USERNAME/NGS-Sequencing-QC-Analysis`

## Option 2: command line

Create an **empty** repository of the same name on GitHub (do NOT initialize it
with a README because the project already has one), then in the extracted folder:

```bash
git init
git add .
git commit -m "Add tested synthetic FASTQ QC Python demonstration"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/NGS-Sequencing-QC-Analysis.git
git push -u origin main
```

Replace `YOUR_USERNAME` with your GitHub username and authenticate if prompted.

## After publishing

- Open README.md on GitHub and confirm both charts render.
- Select **Actions** to review the automated Python tests (may require enabling
  actions in repository settings).
- Check `results/qc_summary.json` and `docs/scientific_notes.md`.
- Link the repository from your Wolf Analytics portfolio only after confirming
  it is visible and you are comfortable with the version.

The package was designed and tested offline; **no changes to your GitHub
account have been made**.
