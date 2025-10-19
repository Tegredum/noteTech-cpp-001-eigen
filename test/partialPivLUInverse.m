% encoding: utf-8
% author: Grok 4 Expert
% matlab version: R2024b

function invA = partialPivLUInverse(A)
	% 输入：方阵 A
	% 输出：A 的逆矩阵 invA
	% 假设 A 可逆且方阵
	
	n = size(A, 1);
	if size(A, 1) ~= size(A, 2)
		error('矩阵必须为方阵');
	end
	
	% 初始化 LU 分解（原地存储于 L 和 U）
	LU = A;  % 复制 A 以进行原地分解
	P = 1:n;  % 初始置换向量（行索引）
	
	% 部分主元 LU 分解
	for k = 1:n-1
		% 寻找当前列 k 中绝对值最大的主元
		[~, pivotIdx] = max(abs(LU(k:n, k)));
		pivotIdx = pivotIdx + k - 1;  % 调整为全局索引
		
		% 如果主元为零，则矩阵奇异（但假设可逆，此处跳过检查）
		
		% 交换行（置换）
		if pivotIdx ~= k
			% 交换 LU 中的行
			tempRow = LU(k, :);
			LU(k, :) = LU(pivotIdx, :);
			LU(pivotIdx, :) = tempRow;
			
			% 交换置换向量 P
			tempP = P(k);
			P(k) = P(pivotIdx);
			P(pivotIdx) = tempP;
		end
		
		% 消除下三角部分（计算 L 和 U）
		for i = k+1:n
			LU(i, k) = LU(i, k) / LU(k, k);  % L 的元素（单位对角）
			for j = k+1:n
				LU(i, j) = LU(i, j) - LU(i, k) * LU(k, j);  % 更新 U
			end
		end
	end
	
	% 现在 LU 包含 L（下三角，非对角部分）和 U（上三角，包括对角）
	
	% 求解 A X = I，即 X = inv(A)
	I = eye(n);  % 单位矩阵
	invA = zeros(n, n);
	
	% 对 I 的每一列求解（使用前向和后向替换）
	for col = 1:n
		b = I(:, col);  % I 的第 col 列
		
		% 应用置换：b = P^{-1} * b（等价于 b(P)）
		bPerm = b(P, :);  % 因为 P 是行置换向量
		
		% 前向替换：求解 L y = bPerm（L 是单位下三角）
		y = zeros(n, 1);
		for i = 1:n
			y(i) = bPerm(i);
			for j = 1:i-1
				y(i) = y(i) - LU(i, j) * y(j);
			end
		end
		
		% 后向替换：求解 U x = y（U 是上三角）
		x = zeros(n, 1);
		for i = n:-1:1
			x(i) = y(i);
			for j = i+1:n
				x(i) = x(i) - LU(i, j) * x(j);
			end
			x(i) = x(i) / LU(i, i);
		end
		
		invA(:, col) = x;
	end
end
