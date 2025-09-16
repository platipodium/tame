<!--
SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH

SPDX-License-Identifier: Apache-2.0
-->

(1)
converts biological unit d-1 into physical FABM/driver unit s-1 for RHS!
physical drivers usually have s-1 as a natural time unit.  So we
need to convert the biological time unit d-1 into s-1.
We agree to make all rate constants and parameters in fabm.yaml in d-1.
Conversion occurs at fabm host interaciton (_ADD_SOURCE_)
We use `secs_per_day` and `days_per_sec` for this conversion.

(2)
use functions only for
 (a) long, but not essential code parts
 (b) if dependency is used more than once

 (3)
