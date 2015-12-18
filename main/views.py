from django.shortcuts import render

def home(request):
    #import subprocess
    #subprocess.call(["hdr_create -s 4 -o /Users/kent/Documents/OnlineHDR/images/hdr.exr /Users/kent/Documents/OnlineHDR/images/image?.jpg"], shell=True)
    return render(request, 'index.html', {})