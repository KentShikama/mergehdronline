host = "http://localhost:8000"
dark_file_path = 'images/kluki_dark.jpg'
normal_file_path = 'images/kluki_normal.jpg'
bright_file_path = 'images/kluki_bright.jpg'

from poster.streaminghttp import register_openers
from poster.encode import multipart_encode
import urllib2, urllib
import json

register_openers()
data, headers = multipart_encode({'dark': open(dark_file_path), 'normal': open(normal_file_path), 'bright': open(bright_file_path), 'stops': 4, 'alpha': 0.1, 'beta': 0.1, 'theta': 0})
request = urllib2.Request(host + "/api/", data, headers)
response = urllib2.urlopen(request).read()
id = json.loads(response)['id']
remote_file_name = "hdr_final_" + str(id) + ".jpg"
urllib.urlretrieve(host + "/media/" + remote_file_name, "sample_out.jpg")