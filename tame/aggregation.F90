! SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
!
! SPDX-License-Identifier: Apache-2.0

! --- phytoplankton aggregation -------------------------------------------------
! If biovolume is primarily determined by the nitrogen content, also for detritus

! aggreg_rate = self%phi_agg * dom%C * (phy%N + det%N)                    ! [d^{-1}]

! call particle_aggregation(self,sens,env)
! aggreg_rate = 0.0_rk
!! dom_dep     = self%agg_doc*dom%C/(1.0_rk+self%agg_doc*dom%C)
!! aggreg_rate = min(self%phi_agg * dom_dep * (phy%N + det%N),0.1_rk)
!         vS * exp(-4*phys_status )                ! [d^{-1}]
! additional mortality due to H2S stress (EC_50 :55 mmol-H2S/m3; Kuester2005  or for Skel cost 3 mmol-H2S/m3, Breteler1991)
!!if (self%BioOxyOn) then
!!   aggreg_rate = aggreg_rate + self%mort_ODU* env%odu
!!endif
