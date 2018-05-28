import tg
from tg import AppConfig
from tg import TGController
from tg import expose
import kajiki


page = kajiki.XMLTemplate(u'''<html>
    <head></head>
    <body>
      <div id="isomor"></div>

      <script py:for="m in g.webassets['bundle.js'].urls()"
              src="$m">
      </script>
      <script>
ReactDOM.render(
    React.createElement(HelloWorld.HelloWorld, { name: "World" }),
    document.getElementById('isomor')
);
      </script>
    </body>
</html>
''', mode='html5')


class RootController(TGController):
    @expose()
    def index(self):
        return page(dict(
            g=tg.app_globals
        )).render()


config = AppConfig(minimal=True, root_controller=RootController())
config.renderers = ['kajiki']
config.serve_static = True
config.paths['static_files'] = 'static'

from webassets.filter import register_filter
from dukpy.webassets import BabelJSX
register_filter(BabelJSX)

import tgext.webassets as wa
wa.plugme(
    config,
    options={
        'babel_modules_loader': 'umd'
    },
    bundles={
        'bundle.js': wa.Bundle(
            'js/react.js',
            'js/react-dom.js',
            wa.Bundle(
                'js/HelloWorld.jsx',
                filters='babeljsx',
            ),
            output='assets/bundle.js'
        )
    }
)

application = config.make_wsgi_app()

from wsgiref.simple_server import make_server
print("Serving on port 8080...")
httpd = make_server('', 8080, application)
httpd.serve_forever()