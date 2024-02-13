"""Scheduler environment module for NotreDameCRC."""

import flow

from .ndcrc_scheduler import SGEScheduler


class NotreDameCRC(flow.environment.StandardEnvironment):
    """Scheduler environment for NotreDameCRC."""

    hostname_pattern = r".*\.crc\.nd\.edu$"
    template = "crc.nd.sh"
    scheduler_type = SGEScheduler
    JOB_ID_SEPARATOR = "-"


# class NotreDameCRCTest(flow.environment.TestEnvironment):
#
#    hostname_pattern = r'.*\.crc\.nd\.edu$'
#    template = 'crc.nd.sh'
