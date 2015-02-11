# """
# Adapted from
#
#   Eric S. Raymond
#
#   Greatly modified by Nigel W. Moriarty
#   April 2003
#
#  PEXPECT LICENSE
#
#      This license is approved by the OSI and FSF as GPL-compatible.
#          http://opensource.org/licenses/isc-license.txt
#
#      Copyright (c) 2012, Noah Spurrier <noah@noah.org>
#      PERMISSION TO USE, COPY, MODIFY, AND/OR DISTRIBUTE THIS SOFTWARE FOR ANY
#      PURPOSE WITH OR WITHOUT FEE IS HEREBY GRANTED, PROVIDED THAT THE ABOVE
#      COPYRIGHT NOTICE AND THIS PERMISSION NOTICE APPEAR IN ALL COPIES.
#      THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
#      WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
#      MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#      ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#      WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#      ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#      OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#  """

import string
import logging
from time import time

import pexpect


log = logging.getLogger('logger.main')
loghost = logging.getLogger('logger.host')


#===============================================================================
#  SSHSession
#===============================================================================
class SSHSession(object):
    def __init__(self, user, host, name="", password=None, verbose=True):
        self.terminated = False
        self.connected = False
        self.user = user
        self.host = host
        self.verbose = verbose
        self.password = password
        self.keys = [
            'authenticity',
            'password:',
            '@@@@@@@@@@@@',
            'Command not found.',
            pexpect.EOF,
        ]
        self.name = name if name != "" else host
        self.history = ["init", ]
        self.timeout = self.init_timeout = 15
        self.connected = self._connect()

    def _connect(self, repetition=5):
        show_time = False
        saveverb = self.verbose
        connected = False
        st = st0 = time()
        log.debug("Connecting to server %s..." % self.name)
        self.verbose = False
        connected = self.is_connected()
        st = time()
        if show_time:
            print st - st0
        if not connected:
            for i in range(repetition):
                self.timeout = self.timeout * 1.1
                self.init_timeout = self.timeout
                connected = self.is_connected()
                st1 = st
                st = time()
                if show_time:
                    print st - st0, st - st1
                if connected:
                    break
        self.verbose = saveverb
        if connected:
            log.info("Connected to server %s" % self.name)
            return True
        else:
            log.info("We cannot connect to the server for more than 30s.")
            log.info("Please, check your network.")
            log.info("Terminated")
            self.terminated = True
            return False

    def is_connected(self):
        out = self.ssh("echo ready", verbose=False)
        error = (out == "ERROR")
        unreachable = ("Network is unreachable" in out)
        ready = (out.strip() == 'ready')
        conn = (not error) and (not unreachable) and ready
        return conn

    def check(self):
        self.connected = self.is_connected()

    def __repr__(self):
        if self.connected:
            return 'SSH2 session host:%(host)s user:%(user)s' % {
                'host': self.host,
                'user': self.user}
        else:
            return ""

    def _exec(self, server, command, verbose=True):
        self.timeout = self.init_timeout  # reinit timeout
        return self.__do_exec(server, command, verbose)

    def __do_exec(self, server, command, verbose):
        #  Execute a command on the remote host.    Return the output.
        com = "%s %s" % (server, command)
        child = pexpect.spawn(
            com.strip(),
            timeout=self.timeout
        )
        if self.verbose and verbose:
            if server != "":
                loghost.debug("%s>%s" % (self.name, command))
            else:
                loghost.debug("localhost>%s" % command)

        try:
            child.expect(self.keys)
            out = child.before
        except:
            out = "ERROR"
            # try again with a slightly increased timeout
            # (limiting the increase to 30s - avoiding to much recursivity!)
            self.timeout = self.timeout * 1.1
            if self.timeout <= 30.0:
                out = self.__do_exec(server, command, verbose)

        self.history.append(out)
        if self.verbose and verbose:
            loghost.debug(out)
        return out

    def ssh(self, command, verbose=True):
        out = self._exec(server="ssh -l %s %s" % (self.user, self.host),
                         command="%s" % command, verbose=verbose)
        return out

    def scp(self, f1, mode, f2):
        """usage: scp(file1, mode=['to','from'], file2)"""
        if mode == 'to':
            return self._exec(server="",
                              command="scp %s %s@%s:%s" % (f1, self.user,
                                                           self.host, f2))
        elif mode == 'from':
            return self._exec(server="",
                              command="scp %s@%s:%s %s" \
                                      % (self.user,
                                         self.host, f2, f1))
        else:
            log.error(
                "%s>scp mode unknown (should be 'to' or 'from')" % self.name)
            log.error(
                "%s>usage: scp(file1, mode=['to','from'], file2)\n" % self.name)

    def exists(self, _file):
        """
        Retrieve file permissions of specified remote file.
        """
        seen = self.ssh("/bin/ls -ld %s" % _file, verbose=False)
        if string.find(seen, "No such file") > -1:
            return False  # File doesn't exist
        else:
            try:
                return seen.split()[0]  # Return permission field of listing.
            except:
                pass


if __name__ == "__main__":

    import sys

    log.addHandler(logging.StreamHandler())
    log.setLevel(logging.INFO)

    s = SSHSession("nmr", "192.168.223.175", verbose=True)
    if not s.connected:
        sys.exit()
    log.info("connected")
    #s.scp("__init__.py", "to", "essai_scp")
    #s.exists("chabazite-0")
    #s.ssh('cd data/siblanixenon/nmr; ls')

    # EOF --------------------------------------------------------------------------