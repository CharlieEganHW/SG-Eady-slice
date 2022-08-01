function l=get_edgeLengths(v)

l=vecnorm(v-circshift(v,-1),2,2);

end
