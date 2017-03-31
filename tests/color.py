
import re
import shlex
from nose.plugins.base import Plugin


def color_arguments(arguments):
    categories = (
        (('l', 'alname', 'alhead', 'allink', 'altail'), '33'),
        (('box', 'pbc', 'x', 'y', 'z'), '32'),
        (('f', ),  '31'),
        (('d', 'dz'), '34'),
    )
    template_from = "(-({}) [a-zA-Z0-9:,-]+)(\\s|$|')"
    template_to = '\\033[{}m\\1\\033[0m\\3'
    for option, color in categories:
        arguments = re.sub(template_from.format('|'.join(option)),
                           template_to.format(color),
                           arguments)
    return arguments


class ColorVerbose(Plugin):
    enabled = True
    env_opt = 'NOSE_COLOR_VERBOSE'
    name = 'color-verbose'
    score = 1600

    def options(self, parser, env):
        """Registers the commandline option, defaulting to enabled.
        """
        parser.add_option(
            "--no-{}".format(self.name), action="store_false",
            default=not env.get(self.env_opt), dest="color_verbose",
            help="Color the tests by category [{}]".format(self.env_opt))

    def configure(self, options, conf):
        """Configure plugin. Plugin is enabled by default.
        """
        self.config = conf
        try:
            self.enabled = options.color_verbose
        except AttributeError:
            self.enabled = False

    def describeTest(iself, test):
        args = test.id()
        return color_arguments(args)

    def report(self, stream):
        stream.write('Color is beautiful\n')
