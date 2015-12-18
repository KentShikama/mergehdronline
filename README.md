# HDR GIMP Plugin

A [HDR merging GIMP plugin](kent_hdr.py) along with an API service that allows you to instantly merge HDRs online for free.

## How to use GIMP plugin

- Install [GIMP](https://www.gimp.org)
- Download the [GIMP plugin (kent_hdr.py)](kent_hdr.py) from this repository and place it in the plug-in directory (e.g., /Users/me/Library/Application\ Support/GIMP/2.8/plug-ins/kent_hdr.py)
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
