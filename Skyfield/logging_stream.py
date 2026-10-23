"""Line-buffered text redirection into an existing logger."""


class StreamToLogger:
    """Collect text writes into complete log lines, flushing partial lines on demand."""

    def __init__(self, logger, log_level):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        self.linebuf += buf
        while '\n' in self.linebuf:
            line, self.linebuf = self.linebuf.split('\n', 1)
            self.logger.log(self.log_level, line.rstrip('\r'))
        return len(buf)

    def flush(self):
        if self.linebuf:
            line, self.linebuf = self.linebuf, ''
            self.logger.log(self.log_level, line.rstrip('\r'))
