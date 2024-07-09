'''
Generate an earth system model configuration from the given configuration file.

$Id: expconfig.py 4103 2015-07-16 07:23:05Z monika $
'''

import os
import re
import StringIO

from itertools import dropwhile

from configobj import ConfigObj, InterpolationError

import feedback

class ExpConfigError(InterpolationError):
    def __init__(self, message, key):
        message = message.rstrip('.!')
        InterpolationError.__init__(self, 
            "{0} while reading key '{1}'".format(message, key))

class ExpConfig(ConfigObj):
    '''Read and store configuration info from input and experiments' library

    Store environment as default for control settings, then add config from files
    '''
    
    #
    # Basic settings
    #

    exp_lib_dir = 'standard_experiments'
    env_lib_dir = 'standard_environments'
    default_name = 'DEFAULT'
    id_name = 'EXP_ID'

    # Class constructor

    def __init__(self, experiment_config_name, extra_dict={}, config_roots=['']):
        '''Read experiment config to get basic settings
        
        TODO: probably nicer if default experiment is given as argument
        '''

        #
        # Helper functions
        #

        def split_jobs(config):
            '''Post-process job definition to allow for shared configs as [[job1, job2]]'''
            if 'jobs' in config:
                sep = re.compile(r'\s*,\s*')
                for subjobs, subconfig in config['jobs'].iteritems():
                    if re.search(sep, subjobs):
                        for subjob in re.split(sep, subjobs):
                            if subjob in config['jobs']:
                                config['jobs'][subjob].merge(subconfig.dict())
                            else:
                                config['jobs'][subjob] = subconfig.dict()
                        del config['jobs'][subjobs]

        def get_config_name(lib_name, base_name):
            '''Cycle through config path until a match is found.
               
               Return simple path otherwise'''
            config_name = os.path.join(lib_name, base_name)
            for config_root in config_roots:
                tentative_name = os.path.join(config_root, config_name)
                if os.path.exists(tentative_name):
                    config_name = tentative_name
                    break
            return config_name

        #
        # Method body
        #

        # Pre-read basic experiment settings

        pre_config = ConfigObj(experiment_config_name, interpolation=False)

        experiment_type = extra_dict.get('EXP_TYPE', pre_config['EXP_TYPE'])
        environment = extra_dict.get('ENVIRONMENT', 
                      pre_config.get('ENVIRONMENT',
                      ExpConfig.default_name))
        # Backwards compatibility ENVIRONMENT -> QUEUE_TYPE
        if environment == ExpConfig.default_name and 'QUEUE_TYPE' in pre_config:
            feedback.warning("found obsolete keyword 'QUEUE_TYPE'; "
                             "should be replaced by 'ENVIRONMENT'")
            environment = pre_config['QUEUE_TYPE']

        pre_config = None

        # Start from empty configuration

        pre_config = ConfigObj(interpolation=False)
        config_versions = []

        # Get default experiment id from file name
        pre_config[ExpConfig.id_name] = os.path.splitext(
            os.path.basename(experiment_config_name)
        )[0]

        # Read Environment

        env_dict = dict(os.environ)
        # Mask literal dollar characters
        for key, value in env_dict.iteritems():
            env_dict[key] = value.replace('$', '$$')
        pre_config.merge({'DEFAULT': env_dict})

        # Read experiment settings from library (default and type specific)

        lib_config_name = get_config_name(ExpConfig.exp_lib_dir,
                                          ExpConfig.default_name+'.config')
        pre_config.merge(ConfigObj(lib_config_name, interpolation=False))
        split_jobs(pre_config)

        lib_config_name = get_config_name(ExpConfig.exp_lib_dir, 
                                          experiment_type+'.config')
        if os.path.exists(lib_config_name):
            pre_config.merge(ConfigObj(lib_config_name, interpolation=False))
            split_jobs(pre_config)
        else:
            feedback.warning("cannot find config for '%s', using default only",
                             experiment_type)

        # Read host environment settings from library

        lib_config_name = get_config_name(ExpConfig.env_lib_dir,
                                          environment+'.config')

        if os.path.exists(lib_config_name):
            pre_config.merge(ConfigObj(lib_config_name, interpolation=False))

        # Re-read config to allow overriding default settings
        # TODO: probably nicer if default experiment is given as argument
        experiment_config = ConfigObj(experiment_config_name,
                                      interpolation=False)
        pre_config.merge(experiment_config)
        split_jobs(pre_config)

        # Add extra dictionary
        pre_config.merge(extra_dict)

        if not 'ENVIRONMENT' in pre_config:
            pre_config['ENVIRONMENT'] = environment

        # Add complete versioning info

        # Re-read merged config with interpolation set.
        # This works around incomprehensible inheritance of interpolation with
        # merge. Make sure that all values are interpolated

        config_lines = StringIO.StringIO()

        pre_config.write(config_lines)
        pre_config = None

        config_lines.seek(0)
        pre_config = ConfigObj(config_lines, interpolation='template')

        # Extract experiment description from initial comment
        # if not set explicitly
        if not pre_config.has_key('EXP_DESCRIPTION'):
            is_empty = lambda s: re.match(r'^\s*$', s)
            rm_comment = lambda s: re.sub(r'^\s*# ?', '', s)       
            pre_config['EXP_DESCRIPTION'] = "\n".join(
                reversed(list(
                    dropwhile(is_empty,
                        reversed(list(
                            dropwhile(is_empty,
                                map(rm_comment,
                                    experiment_config.initial_comment)
                            )
                        )) 
                    )
                ))
            )

        def read_value(value):
            if os.path.exists(value):
                stream = open(value)
                result = stream.read().strip()
                stream.close()
            else:
                result = ''
            return result

        def sec2time(seconds):
            '''Create time string (HH:MM:SS) from second of day'''
            seconds = int(seconds)
            if seconds >= 86400:
                raise ValueError("invalid second of day '{0}'".format(seconds))
            minutes, s = divmod(seconds, 60)
            h, m = divmod(minutes, 60)
            return "{0:02}:{1:02}:{2:02}".format(h, m, s)

        def split_date(value):
            '''Re-format datetime string to list for use in namelists'''
            match = re.match(r'^0*(\d+)-0*(\d+)-0*(\d+)'
                             r'([T ]0*(\d+)(:0*(\d+)(:0*(\d+))?)?)?$', value)
            if match:
                numbers = match.groups('0')
            else:
                raise ValueError("invalid date/time '{0}'".format(value))
            return [numbers[i] for i in [0,1,2,4,6,8]]

        def eval_value(value):
            '''
                Evaluate key as python expression,
                return as string or sequence of strings.
            '''
            result = eval(value)
            if isinstance(result, (list, tuple)):
                result = map(str, result)
            else:
                result = str(result)
            return result

        def eval_value_string(value):
            '''
                Evaluate key as python expression,
                return as string or sequence of strings.
            '''
            result = eval_value(value)
            if isinstance(result, (list, tuple)):
                result = ", ".join(result)
            return result

        def eval_expression(value):
            '''
                Check if value is a supported expression.
                If so, evaluate and return result, otherwise just pass through.
            '''
            match = re.match(r'^eval\((.*)\)$', value, re.S)
            if match:
                return eval_value(match.group(1))

            match = re.match(r'^evals\((.*)\)$', value, re.S)
            if match:
                return eval_value_string(match.group(1))

            match = re.match(r'^split_date\((.*)\)$', value, re.S)
            if match:
                return split_date(match.group(1))

            match = re.match(r'^sec2time\((.*)\)$', value, re.S)
            if match:
                return sec2time(match.group(1))

            match = re.match(r'^read\((.*)\)$', value, re.S)
            if match:
                return read_value(match.group(1))

            return value

        # Interpolate and evaluate keys if they are an expression
        def eval_key(section, key):
            try:
                value = section[key]
                if isinstance(value, (list, tuple)):
                    value = map(eval_expression, value)
                elif isinstance(value, basestring):
                    value = eval_expression(value)                    
            except (InterpolationError, ValueError) as error:
                raise ExpConfigError(error.message, key)
            section[key] = value

        pre_config.walk(eval_key)

        # Re-read final config without interpolation.
        # This allows copying data without evaluation of version keywords.

        config_lines.seek(0)
        config_lines.truncate()

        pre_config.write(config_lines)
        pre_config = None

        config_lines.seek(0)
        ConfigObj.__init__(self, config_lines, interpolation=False)
        
        self.experiment_id = self[ExpConfig.id_name]
        self.experiment_kind = re.sub(r'-\w+$', '', experiment_type)

