from api.views import FileUploadView
from django.conf.urls import url

urlpatterns = [
    url(r'^$', FileUploadView.as_view(), name="upload"),
]