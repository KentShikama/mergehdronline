from django.core.files.storage import FileSystemStorage
import os
import random
from rest_framework.parsers import FileUploadParser
from rest_framework.views import APIView
from rest_framework.response import Response

class FileUploadView(APIView):
    parser_classes = (FileUploadParser,)

    def cleanup(self, bright_path, dark_path, exr_path, fs, normal_path, tonemapped_hdr_path):
        fs.delete(dark_path)
        fs.delete(normal_path)
        fs.delete(bright_path)
        fs.delete(exr_path)
        fs.delete(tonemapped_hdr_path)

    def merge_hdr_as_subprocess(self, MEDIA_PATH, bright_path, dark_path, exr_path, final_jpg_path, normal_path,
                                tonemapped_hdr_path, stops, alpha, beta, theta):
        import subprocess
        hdr_create_command = "hdr_create -s " + stops + " -o " + exr_path + " " + os.path.join(MEDIA_PATH,
                dark_path) + " " + os.path.join(MEDIA_PATH, normal_path) + " " + os.path.join(MEDIA_PATH, bright_path)
        hdr_tonemap_command = "hdr_squish -a " + alpha + " -b " + beta + " -t " + theta +  "  -i " + exr_path + " -o " + tonemapped_hdr_path
        hdr_convert_command = "hdr_convert -i " + tonemapped_hdr_path + " -o " + final_jpg_path
        subprocess.call([hdr_create_command], shell=True)
        subprocess.call([hdr_tonemap_command], shell=True)
        subprocess.call([hdr_convert_command], shell=True)

    def save_hdr_photos(self, request):
        dark = request.FILES['dark']
        normal = request.FILES['normal']
        bright = request.FILES['bright']
        id = str(random.randint(1, 999999999))
        fs = FileSystemStorage()
        dark_path = fs.save('dark_' + id + '.jpg', dark)
        normal_path = fs.save('normal_' + id + '.jpg', normal)
        bright_path = fs.save('bright_' + id + '.jpg', bright)
        return bright_path, dark_path, fs, id, normal_path

    def create_output_file_paths(self, id):
        BASE_PATH = os.path.dirname(os.path.dirname(__file__))
        MEDIA_PATH = os.path.join(BASE_PATH, "HDR", "media")
        exr_name = "hdr_" + id + ".exr"
        exr_path = os.path.join(MEDIA_PATH, exr_name)
        tonemapped_hdr_name = "hdr_tonemapped_" + id + ".exr"
        tonemapped_hdr_path = os.path.join(MEDIA_PATH, tonemapped_hdr_name)
        final_jpg_name = "hdr_final_" + id + ".jpg"
        final_jpg_path = os.path.join(MEDIA_PATH, final_jpg_name)
        return MEDIA_PATH, exr_path, final_jpg_path, tonemapped_hdr_path

    def post(self, request, format=None):
        stops = request.data['stops']
        alpha = request.data['alpha']
        beta = request.data['beta']
        theta = request.data['theta']
        bright_path, dark_path, fs, id, normal_path = self.save_hdr_photos(request)
        MEDIA_PATH, exr_path, final_jpg_path, tonemapped_hdr_path = self.create_output_file_paths(id)
        self.merge_hdr_as_subprocess(MEDIA_PATH, bright_path, dark_path, exr_path, final_jpg_path, normal_path, tonemapped_hdr_path, stops, alpha, beta, theta)
        self.cleanup(bright_path, dark_path, exr_path, fs, normal_path, tonemapped_hdr_path)
        data = {'id': id }
        return Response(status=200, data=data)