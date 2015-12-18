from django import template
register = template.Library()

@register.simple_tag
def active_nav(request, word):
    if not word and request.path == "/":
        return "active"
    elif word and word in request.path:
        return "active"
    else:
        return ""