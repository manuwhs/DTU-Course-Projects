from os import path
from webassets import Environment, Bundle

try:
    from urllib.parse import urlparse
except ImportError:
    from urlparse import urlparse

from tg.configuration.utils import coerce_config
from tg.support.converters import asbool, asint, aslist

import logging
log = logging.getLogger('tgext.webassets')


class SetupExtension(object):
    def __init__(self, configurator, options):
        self.configurator = configurator
        self.options = options

    def __call__(self):
        from tg import config
        # Save as a class attribute so it is available between multiple app instances
        config['tg.app_globals'].__class__.webassets = self.make_webassets_env_from_config(self.options,
                                                                                           config)
    def string_or_bool(self, value):
        try:
            lowervalue = value.lower()
        except AttributeError:
            return asbool(value)

        if lowervalue in ('false', 'no', 'off', 'n', 'f', '0',
                          'true', 'yes', 'on', 'y', 't', '1'):
            return asbool(value)

        return value

    def make_webassets_env_from_config(self, options, config):
        settings = {}
        settings.update(options)
        settings.update(coerce_config(config, 'webassets.', {
            'debug': asbool,
            'auto_build': asbool,
            'manifest': self.string_or_bool,
            'url_expire': asbool,
            'cache': self.string_or_bool,
            'load_path': aslist
        }))

        static_files_path = config['paths']['static_files']
        asset_dir = settings.pop('base_dir', static_files_path)

        asset_url = settings.pop('base_url', '/')
        if not asset_url.startswith('/'):
            if urlparse(asset_url).scheme == '':
                asset_url = '/' + asset_url

        assets_env = Environment(asset_dir, asset_url, **settings)

        bundles = self.options.get('bundles', {})
        for name, value in bundles.items():
            if isinstance(value, Bundle):
                assets_env.register(name, value)
            else:
                # If it's not a bundle consider it a loader
                assets_env.register(value.load_bundles())

        return assets_env
