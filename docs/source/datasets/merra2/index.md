# MERRA2

This page explains how you can setup access MERRA2 data from NASA's [GES DISC](https://disc.gsfc.nasa.gov/).

## Creating an Earthdata Login Profile and Approving the GES DISC App

To download MERRA2 data, you'll need to create an [Earthdata Login](https://urs.earthdata.nasa.gov/users/new).

After creating an account, you'll need to connect your account to GES DISC to access data.

The following procedure is sourced from [GES DISC's documentation](https://disc.gsfc.nasa.gov/earthdata-login):
1. Login to [https://urs.earthdata.nasa.gov/](https://urs.earthdata.nasa.gov/).

2. Click on the 'Applications' tab in the top menu bar, and in the pop up menu, select the 'Authorized Apps' tab.

3. Scroll to bottom of page to find the "Approve More Applications" button. Click it to open the search bar.

4. Search for and select the **NASA GESDISC DATA ARCHIVE**.  Once you've clicked approve, agree to the EULA.

5. The data source has now been connected to your account.

## Configure API Crendentials

To download MERRA2 data via **geodata** you'll need to install the API credentials locally.  

### macOS/Linux

For MacOS, create a file in your root directory (`cd ~`) called `.netrc` by opening a terminal window and running `touch ~/.netrc` to create the file and then execute the following command:

```bash
echo "machine urs.earthdata.nasa.gov login [login] password [password] " >> .netrc
```

where `[login]` is your Earthdata user name and `[password]` is your Earthdata Login password.

### Windows

For Windows, open Notepad and enter the following line in a new document, making sure to substitute `<uid>` and `<password>`for your Earthdata login credentials:
 
`machine urs.earthdata.nasa.gov login <uid> password <password>`

Save the file to `C:\Users\<username>\.netrc`

## What' next?

We provide a few more tutorials on how to download and utilize the MERRA2 dataset with Geodata. You can find them below.

```{toctree}
:maxdepth: 2
:glob:
*
```
