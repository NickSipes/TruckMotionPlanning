function semitruck_steer(pub_l,pub_r,msg_l,msg_r,pos)
    msg_l.Data = pos;
    msg_r.Data = pos;
    send(pub_l,msg_l);
    send(pub_r,msg_r);
end