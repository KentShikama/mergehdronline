# No longer maintained

# HDR GIMP Plugin

A [HDR merging GIMP plugin](kent_hdr.py) along with an online HDR merging API service.

Please take a look at [the recording demo](live_demo_of_before_and_after.mov) which first shows the feca's plugin and then shows mine. You can download the entire repository by clicking the Download ZIP button on Github.

## How to use GIMP plugin

- Install [GIMP](https://www.gimp.org)
- Download the [GIMP plugin (kent_hdr.py)](kent_hdr.py) from this repository and place it in the plug-in directory (e.g., /Users/me/Library/Application\ Support/GIMP/2.8/plug-ins/kent_hdr.py)
  - Please see the [GIMP wiki for generic instructions on how to install GIMP plugins](https://en.wikibooks.org/wiki/GIMP/Installing_Plugins)
- Load the exposure bracketed photos into GIMP (Normal, -EV, +EV) as 3 layers
- From the menu: Python-Fu -> Merge with Gradient Based HDR
- Adjust the parameters if you wish (the defaults should most be fine)

![Options](images/plugin_options_screen.png?raw=true)

- Click OK

## Motivation 

As of Dec 18, 2015, GIMP doesn't have a prepackaged HDR merging filter. The GIMP plugin repository has a number of HDR plugins but none of them implement a tone mapping operator that preserves local contrasts. I wanted to create a GIMP plugin that would allow users to be able to experience the power of the latest HDR innovations. After a quick search on Google Scholar, the [Gradient Domain High Dynamic Range Compression algorithm by Raanan Fattal, Dani Lischinski and Michael Werman](http://ahtariev.ru/OLD/content/hdr_art/hdrc.pdf) seemed to be the latest HDR tone mapping algorithm with subjectively the best results, and thus this algorithm is what I used in the GIMP plugin and API.

## What is HDR merging? And what is the significance of this HDR GIMP plugin?

Normal photographs can only capture a certain span of the wavelengths that the human eye can capture. Hence, if you take a photo with a lower exposure, you would be able to see details in the brighter regions but the details in the darker regions would become all black, and vice versa. The general idea of HDR merging is that if you take multiple photographs of the same spatial scene and merge them together, you get a photograph that captures both the details in the brighter regions and the details in the darker regions.

There are a couple ways that HDR merging can be implemented. Assume we are working with 3 images: normal exposure, -EV (a darker image), and +EV (a brighter image). An example of three such images is shown below (from http://www.easyhdr.com/examples/):

![Dark](images/kluki_dark.jpg?raw=true)

-EV (Dark exposure) Photograph by Bartłomiej Okonek

![Normal](images/kluki_normal.jpg?raw=true)

Normal exposure Photograph by Bartłomiej Okonek

![Bright](images/kluki_bright.jpg?raw=true)

+EV (Bright exposure) Photograph by Bartłomiej Okonek

One way HDR merging method is to simply lay the 3 images on top of each other and mask out dark regions from the -EV and normal image layers and mask out bright regions from the +EV and normal image layers. And then optionally at the end, scale the image intensities so that it covers the full range of your display monitor. [Feca's HDR (High Dynamic Range) with tone mapping plugin](http://registry.gimp.org/node/24310), which was the inspiration for this plugin, essentially performs this method of HDR merging. Here is the result of using Feca's plugin:

![Feca HDR](images/kluki_using_feca_hdr.jpg?raw=true)

Resulting HDR merged image using Feca's plugin

While this does indeed allow you to preserve details as you can see above, it often causes unwanted halos. This is especially apparent in the sky and the edge of the roof.

The HDR merging method used in this GIMP plugin and API, as mentioned above, instead creates a HDR radiance map which reflects the "true" radiance of the scene by combining the 3 images. Once this HDR radiance map is reconstructed, we scale down the dynamic range of the radiance map (this operation is known as tone mapping). This is needed because as most monitors do not have the capability to display such a dynamic range. While we can simply squish down the dynamic range to fit within the display range of the monitor, this type of "global tone mapping operator" causes local contrasts to disappear. Thus this plugin instead uses an algorithm called [Gradient Domain High Dynamic Range Compression algorithm which is by Raanan Fattal, Dani Lischinski and Michael Werman](http://ahtariev.ru/OLD/content/hdr_art/hdrc.pdf) that squishes larger gradients in the HDR radiance map more than the smaller gradients. This has the effect of preserving local contrasts while still allowing the resulting image to fit within the display range of the monitor. Here is the result of using this plugin (with the default parameters):

![Kent HDR](images/kluki_using_kent_plugin_hdr.jpg?raw=true)

Resulting HDR merged image using this plugin

As you can see, we no longer have the halos in the sky or at the edge of the roof. If you compare this image to the original "normal" exposure photograph, we can see how there is much more details in the sky. Note I made the plugin automatically apply a white balance correction and I have also manually adjusted the brightness and contrast a bit.

## API

The plugin uses an API that I created and am hosting at [Merge HDR Online](http://mergehdronline.com). There are two API endpoints: one for uploading the images you want to merge and one for downloading the resulting JPG.

### API Merge Uploading Endpoint

Make sure to remember the returned ID number as you will need it to retrieve the resulting merged image.

Endpoint (POST): http://mergehdronline.com/api/

Params:

- 'dark': File
- 'normal': File
- 'bright': File
- 'stops': Int (1-8) - Number of stops in your image
- 'alpha': Float (0-0.5) - Strength of preservation of local gradients
- 'beta': Float (0-0.5) - Strength of tone mapping
- 'theta': Float (0-1) - Large scale contrast smoothing

Returns:

- 'id': ID number

### API Resulting JPG Downloading

Endpoint (GET): http://mergehdronline.com/media/hdr_final_{your_id_from_above}.jpg

### Sample API Request/Response Script

You can find [a sample python script in this repository](sample_api_request_script.py) that uploads three images and saves the output to sample_output.jpg. Note that you will need to install poster, which can be done through `pip install poster`.

The API works by wrapping a very minorly tweaked version of the C++ implementation of the Gradient Based HDR compression algorithm by [Andrew Cotter](http://ttic.uchicago.edu/~cotter/projects/hdr_tools/). The web requests are handled by Django. The rationale for using an API is it allows the user to simply install a GIMP plugin in the standard fashion without having to configure and install C scripts or have a performance hit by having the entire HDR merging process executed in python. In decent bandwidth areas, the server should much faster at processing than a local python script even when including time that it takes to POST and GET the files. Note you can also set up the API on your localhost as instructed below in the "Development Set Up" section.

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
