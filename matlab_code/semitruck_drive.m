function semitruck_drive(pub_l,pub_r,msg_l,msg_r,vel)
    msg_l.Data = vel;
    msg_r.Data = vel;
    send(pub_l,msg_l)
    send(pub_r,msg_r)
end