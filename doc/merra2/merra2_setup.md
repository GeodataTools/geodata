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

5. The data source has now been connected to your account.

## Enabling data downloads via the Geodata package

To download MERRA2 data via **geodata** you'll need to install the API credentials locally.  

 * For MacOS, create a file in your root directory (`cd ~`) called `.netrc` by opening a terminal window and running `touch ~/.netrc` to create the file and then:

```
echo "machine urs.earthdata.nasa.gov [login] [password] " >> .netrc
```

 where `[login]` is your Earthdata user name and `[password]` is your Earthdata Login password.

 Alternatively, if you are using VSCode, you can simply open your root folder in a new window, save a new file as `.netrc`, and copy paste in the following:

```
machine urs.earthdata.nasa.gov [login] [password] "
```

* For Windows, open Notepad and enter the following line in a new document, making sure to substitute `<uid>` and `<password>`for your Earthdata login credentials:
 
 `machine urs.earthdata.nasa.gov login <uid> password <password>`

Save the file to `C:\Users\<username>\.netrc`


 Finally, you will also need to create a file called `.urs_cookies` to persist cookies across multiple download calls.  You can create this on Mac OS by opening a terminal window and running:

 ```
 touch ~/.urs_cookies
 ```

in your root directory.

For Windows, open a Powershell window, navigate to ``C:\Users\<username>\` and run the following:

```
echo $null >> .urs_cookies
```


To confirm you have set up MERRA2 data access correctly, see: [MERRA2 Download](https://github.com/east-winds/geodata/blob/master/doc/merra2/merra2_download.md))
