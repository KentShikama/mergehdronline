# HDR GIMP Plugin

A [HDR merging GIMP plugin](kent_hdr.py) along with an API that allows you to instantly merge HDRs online for free.

## How to use GIMP plugin

- Install [GIMP](https://www.gimp.org)
- Download the [GIMP plugin (kent_hdr.py)](kent_hdr.py) from this repository and place it in the plug-in directory (e.g., /Users/me/Library/Application\ Support/GIMP/2.8/plug-ins/kent_hdr.py)
- Load the exposure bracketed photos into GIMP (Normal, -EV, +EV) as 3 layers
- From the menu: Python-Fu -> Merge with Gradient Based HDR
- Adjust the parameters if you wish (the defaults should most be fine)
- Click OK

## Motivation

As of Dec 18, 2015, GIMP doesn't have a prepackaged HDR filter. The GIMP plugin repository has a number of HDR plugins but none of them implement a tone mapping operator that preserves local contrasts. I wanted to create a GIMP plugin that would allow users to be able to experience the power of the latest HDR innovations. After a quick search on Google Scholar, the [Gradient Domain High Dynamic Range Compression algorithm by Raanan Fattal, Dani Lischinski and Michael Werman](http://ahtariev.ru/OLD/content/hdr_art/hdrc.pdf) seemed to be the latest HDR tone mapping algorithm with subjectively the best results, and thus this algorithm is what I used in the GIMP plugin and API.

The gimp plugin was inspired by [feca's HDR (High Dynamic Range) with tone mapping plugin](http://registry.gimp.org/node/24310), although none of the code from feca's plugin remains. The API works by wrapping a very minorly tweaked version of the C++ implementation of the Gradient Based HDR compression algorithm by [Andrew Cotter](http://ttic.uchicago.edu/~cotter/projects/hdr_tools/). 

## Development Set Up

### "HDR Tools" (By Andrew Cotter) Set Up 

See ["HDR Tools" site](http://ttic.uchicago.edu/~cotter/projects/hdr_tools/)

- cd c++/
- Install popt, LibTIFF, OpenEXR and ImageMagick
- Tweak the Makefile to your platform
- make
- make install

### Django Set Up

- Install [Python 3.4+](https://www.python.org/downloads/)
- Copy `HDR/settings.py.example` to `HDR/settings.py` and edit as necessary
- Install all dependencies `pip3 install -r requirements.txt`
- Run the server `python3 manage.py runserver`
- Check out localhost:8000 on your browser
