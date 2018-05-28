from webassets import Bundle

import logging
log = logging.getLogger('tgext.webassets')


def plugme(configurator, options=None, **kwargs):
    if options is None:
        options = {}
    options.update(kwargs)

    log.info('Setting up tgext.webassets extension...')
    from tg.configuration import milestones
    from .setup_extension import SetupExtension
    milestones.environment_loaded.register(SetupExtension(configurator, options))

    return dict(appid='tgext.webassets')


