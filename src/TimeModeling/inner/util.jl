function rotate_Nx2_array(arr, ϕ)
  return arr * [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
end

function rotate_2xN_array(arr, ϕ)
  return [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)] * arr
end
