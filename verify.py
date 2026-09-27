#!/usr/bin/env python3
"""One command: the finite-field checks (check.py), then the direct GR(4, m) lift checks (gr_lift.py).

    python3 verify.py          # about a minute
    python3 verify.py --all    # also evaluates the m = 9 lift on all 262,144 elements of GR(4, 9)

Exits non-zero if any check fails.
"""
import sys

import check
import gr_lift

if __name__ == '__main__':
    print('== finite-field checks (check.py)')
    status = check.main()
    print('\n== direct lift checks on GR(4, m) (gr_lift.py)')
    status |= gr_lift.main()
    print('\nVERIFIED' if status == 0 else '\nFAILED')
    sys.exit(status)
