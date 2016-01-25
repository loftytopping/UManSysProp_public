# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright 2014 Dave Jones <dave@waveform.org.uk>.
#
# This file is part of umansysprop.
#
# umansysprop is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 2 of the License, or (at your option) any later
# version.
#
# umansysprop is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# umansysprop.  If not, see <http://www.gnu.org/licenses/>.

# Python 3 compatibility
from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')

import os
import sys
import binascii
import zipfile
import time

try:
    import zlib # We may need its compression method
    crc32 = zlib.crc32
except ImportError:
    zlib = None
    crc32 = binascii.crc32


ZIP_STORED = zipfile.ZIP_STORED
ZIP_DEFLATED = zipfile.ZIP_DEFLATED


class ZipFile(zipfile.ZipFile):
    def write(self, filename_or_obj, arcname=None, compress_type=None,
            mode=0o100664, modified=None):
        """Put the bytes from filename_or_obj into the archive under the name
        arcname.

        Extended version of write which permits file-like objects to be added
        to zipfiles. If a file-like object is passed as filename_or_obj, the
        mode parameter dictates the mode of the file in the zipfile (defaults
        to regular file, user read-write, group read-write, other read), and
        the modified parameter specifies the file modification time in seconds
        since the UNIX epoch (defaults to now). All other parameters are as in
        the original method.
        """
        if isinstance(filename_or_obj, (str, bytes)):
            return super(ZipFile, self).write(filename_or_obj, arcname, compress_type)
        fp = filename_or_obj
        if modified is None:
            modified = time.time()
        modified = time.localtime(modified)
        if arcname is None:
            if not hasattr(fp, 'name'):
                raise ValueError(
                    'arcname must be specified when filename_or_obj has no '
                    'name attribute')
            arcname = fp.name
        arcname = os.path.normpath(os.path.splitdrive(arcname)[1])
        while arcname[0] in (os.sep, os.altsep):
            arcname = arcname[1:]
        zinfo = zipfile.ZipInfo(arcname, modified[0:6])
        if compress_type is None:
            zinfo.compress_type = self.compression
        else:
            zinfo.compress_type = compress_type
        zinfo.external_attr = (mode & 0xFFFF) << 16
        zinfo.flag_bits = 0x00
        zinfo.CRC = CRC = 0
        zinfo.file_size = 0
        zinfo.compress_size = 0
        zinfo.header_offset = self.fp.tell()

        self._writecheck(zinfo)
        self._didModify = True

        # We've no idea what the file-size is yet, but we need to write the
        # header here. If ZIP64 headers are permitted in the file we'll use
        # them regardless of whether we need them or not
        if sys.version_info >= (2, 7, 4):
            self.fp.write(zinfo.FileHeader(self._allowZip64))
        else:
            self.fp.write(zinfo.FileHeader())
        if zinfo.compress_type == ZIP_DEFLATED:
            cmpr = zlib.compressobj(zlib.Z_DEFAULT_COMPRESSION,
                 zlib.DEFLATED, -15)
        else:
            cmpr = None
        while True:
            buf = fp.read(1024 * 8)
            if not buf:
                break
            zinfo.file_size += len(buf)
            zinfo.CRC = crc32(buf, zinfo.CRC) & 0xFFFFFFFF
            if cmpr:
                buf = cmpr.compress(buf)
                zinfo.compress_size += len(buf)
            self.fp.write(buf)
        if cmpr:
            buf = cmpr.flush()
            zinfo.compress_size += len(buf)
            self.fp.write(buf)
        else:
            zinfo.compress_size = zinfo.file_size
        # Seek backwards and re-write file header (which will now include
        # correct CRC and file sizes)
        position = self.fp.tell()       # Preserve current position in file
        self.fp.seek(zinfo.header_offset, 0)
        if sys.version_info >= (2, 7, 4):
            self.fp.write(zinfo.FileHeader(self._allowZip64))
        else:
            self.fp.write(zinfo.FileHeader())
        self.fp.seek(position, 0)
        self.filelist.append(zinfo)
        self.NameToInfo[zinfo.filename] = zinfo

