# Merge HDRs Online

A [GIMP plugin](kent_hdr.py) along with an API that allows you to instantly merge HDRs online for free.

The C++ implementation of the Gradient Based HDR compression algorithm was created by [Andrew Cotter](http://ttic.uchicago.edu/~cotter/projects/hdr_tools/). The gimp plugin was inspired by [feca's HDR (High Dynamic Range) with tone mapping plugin](http://registry.gimp.org/node/24310).

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