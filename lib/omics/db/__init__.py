"""
This is the omics database module

This is 'stand-alone' Django ORM and db backend.
"""
config = {
    'INSTALLED_APPS': ['omics.db.apps.OmicsDBConfig'],
    'DATABASES': {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': None,
        },
    },
}
