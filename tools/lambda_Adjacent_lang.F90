                      ! if(dlambda.ge.1.0 .and. dlambda.lt.sqrt(2.0) ) then
                      !   real_neighbour = real_neighbour  + 1.0
                      !   neighbours(ip) = neighbours(ip) + 1.0
                      !
                      ! else if (dlambda.ge.sqrt(2.0) .and. dlambda .lt. sqrt(3.0)) then
                      !   ! se esssa situacao verificar-se, precisamos verificar se p e q
                      !   ! possuem algum ponto 6-adjacent
                      !
                      !   ! podemos utilizar um unico loop para verificar os vizinhos de p e q
                      !   adj6(1) = lxyz_inv(lxyz(ip,1),lxyz(ip,2),lxyz(ip,3)+1)
                      !   adj6(2) = lxyz_inv(lxyz(ip,1),lxyz(ip,2),lxyz(ip,3)-1)
                      !   adj6(3) = lxyz_inv(lxyz(ip,1),lxyz(ip,2)+1,lxyz(ip,3))
                      !   adj6(4) = lxyz_inv(lxyz(ip,1),lxyz(ip,2)-1,lxyz(ip,3))
                      !   adj6(5) = lxyz_inv(lxyz(ip,1)+1,lxyz(ip,2),lxyz(ip,3))
                      !   adj6(6) = lxyz_inv(lxyz(ip,1)-1,lxyz(ip,2),lxyz(ip,3))
                      !
                      !   neighbour_tf = .true.
                      !
                      !   do ia = 1, 6
                      !     if(phis(adj6(ia)).gt.0) then
                      !       adjacent6 = anint(real(sum((lxyz(adj6(ia),1:3)-lxyz(ip2,1:3))**2)))
                      !       if(adjacent6.eq.1.00)  then
                      !         neighbour_tf = .false.
                      !       end if
                      !
                      !     end if
                      !   end do
                      !
                      !   if(neighbour_tf) then
                      !     real_neighbour  = real_neighbour  + 1
                      !     neighbours(ip) = neighbours(ip) + 1.0
                      !   end if
                      !
                      !
                      ! else if(dlambda.ge.sqrt(3.0)) then
                      !
                      !   adj18(1) = lxyz_inv(lxyz(ip,1)+1,lxyz(ip,2)+1,lxyz(ip,3)+0)
                      !   adj18(2) = lxyz_inv(lxyz(ip,1)+1,lxyz(ip,2)-1,lxyz(ip,3)+0)
                      !   adj18(3) = lxyz_inv(lxyz(ip,1)-1,lxyz(ip,2)+1,lxyz(ip,3)+0)
                      !   adj18(4) = lxyz_inv(lxyz(ip,1)-1,lxyz(ip,2)-1,lxyz(ip,3)+0)
                      !   adj18(5) = lxyz_inv(lxyz(ip,1)+1,lxyz(ip,2)+0,lxyz(ip,3)-1)
                      !   adj18(6) = lxyz_inv(lxyz(ip,1)+0,lxyz(ip,2)+1,lxyz(ip,3)-1)
                      !   adj18(7) = lxyz_inv(lxyz(ip,1)-1,lxyz(ip,2)+0,lxyz(ip,3)-1)
                      !   adj18(8) = lxyz_inv(lxyz(ip,1)+0,lxyz(ip,2)-1,lxyz(ip,3)-1)
                      !   adj18(9) = lxyz_inv(lxyz(ip,1)+1,lxyz(ip,2)+0,lxyz(ip,3)+1)
                      !   adj18(10) = lxyz_inv(lxyz(ip,1)+0,lxyz(ip,2)+1,lxyz(ip,3)+1)
                      !   adj18(11) = lxyz_inv(lxyz(ip,1)-1,lxyz(ip,2)+0,lxyz(ip,3)+1)
                      !   adj18(12) = lxyz_inv(lxyz(ip,1)+0,lxyz(ip,2)-1,lxyz(ip,3)+1)
                      !
                      !   neighbour_tf = .true.
                      !
                      !   do ia = 1, 18
                      !     if(phis(adj18(ia)).gt.0) then
                      !       adjacent18 = anint(real(sum((lxyz(adj18(ia),1:3)-lxyz(ip2,1:3))**2)))
                      !       if(adjacent18.eq.2.00)  then
                      !         neighbour_tf = .false.
                      !       end if
                      !     end if
                      !   end do
                      !
                      !   if(neighbour_tf) then
                      !     real_neighbour  = real_neighbour  + 1
                      !     neighbours(ip) = neighbours(ip) + 1.0
                      !   end if
                      !
                      ! end if
