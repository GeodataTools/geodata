# MERRA2

How to access MERRA2 data from NASA's [GES DISC](https://disc.gsfc.nasa.gov/).

## Creating an Earthdata Login Profile and Approving the GES DISC App

To download MERRA2 data, you'll need to create an [Earthdata Login](https://urs.earthdata.nasa.gov/users/new).

After creating an account, you'll need to connect your account to GES DISC to access data.

The following procedure is sourced from: https://disc.gsfc.nasa.gov/earthdata-login
1. Login to https://urs.earthdata.nasa.gov/. 

2. Click on the 'Applications' tab in the top menu bar, and in the pop up menu, select the 'Authorized Apps' tab.

3. Scroll to bottom of page to find the "Approve More Applications" button. Click it to open the search bar.

4. Search for and select the "NASA GESDISC DATA ARCHIVE".  Once you've clicked approve, agree to the ELUA.

5. The data source has now been connected to your account and you can begin downloading MERRA2 data through the *geodata* package.

To confirm access, see the section on Creating Cutouts with MERRA2 data (link here).