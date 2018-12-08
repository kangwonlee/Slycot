import configparser
import os


def get_file_path():
    return os.path.split(__file__)[0]


def notice(s):
    print("%s %s" % (__file__, s))


def main():

    changes = 0

    if not os.path.exists(os.path.join(get_file_path(), 'control')):
        notice('unable to find control')
    elif not os.path.isdir(os.path.join(get_file_path(), 'control')):
        notice('control not a path')
    elif not os.path.exists(os.path.join(get_file_path(), 'control', 'setup.cfg')):
        notice('unable to find control/setup.cfg')
    else:
        setup_cfg = configparser.ConfigParser()
        setup_cfg.read(os.path.join(get_file_path(), 'control', 'setup.cfg'))
        if setup_cfg.has_section('nosetests'):
            notice('already has section nosetets')
        else:
            notice('does not have section nosetets yet')
            setup_cfg.add_section('nosetests')
            changes += 1

        if not setup_cfg.has_option('nosetests', 'nocapture'):
            notice('does not have nosetets.nocapture yet')
            setup_cfg.set('nosetests', 'nocapture', '1')
            changes += 1
        else:
            existing_nocapture = setup_cfg.get('nosetests', 'nocapture')
            notice('already has nosetets.nocapture %s' % existing_nocapture)
            if '0' == existing_nocapture:
                notice('setting nosetets.nocapture to 1')
                setup_cfg.set('nosetests', 'nocapture', '1')
                changes += 1

        if changes:
            notice('writing setup.cfg file')
            with open (os.path.join(get_file_path(), 'control', 'setup.cfg'), 'w') as f:
                setup_cfg.write(f)


if "__main__" == __name__:
    main()
