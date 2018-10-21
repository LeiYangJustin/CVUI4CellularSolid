#include "skeleton_extractor.h"

CSkeletonExtractor::CSkeletonExtractor(CImgData img_data):has_void_skeleton_(false)
{
	img_data_ = img_data;
}

CSkeletonExtractor::~CSkeletonExtractor()
{
}


std::vector<cv::Point> CSkeletonExtractor::getSolidSkeleton()
{
	extract_morphological_skeleton(true);
	return solid_skeleton_;
}

std::vector<cv::Point> CSkeletonExtractor::getVoidSkeleton()
{
	has_void_skeleton_ = true;
	extract_morphological_skeleton(false);
	return void_skeleton_;
}

void CSkeletonExtractor::getVoidSkeletonSamples(std::vector<cv::Point>& X)
{
	if (!has_void_skeleton_)
		CSkeletonExtractor::getVoidSkeleton();

	// get void skeleton samples as the local extrema

	cv::Mat void_skel_img;
	convert_skeleton_to_img(void_skeleton_, void_skel_img);
	flood_fill_sampling(void_skel_img, X);
}

void CSkeletonExtractor::extract_morphological_skeleton(bool is_solid)
{
	// 
	cv::Mat skelImg, img;
	if (is_solid) {
		img_data_.GetSolidImg().copyTo(img);
	}
	else {
		img_data_.GetVoidImg().copyTo(img);
	}

	// Skeletonization simple method
	if (false) {
		cv::threshold(img, img, 127, 255, cv::THRESH_BINARY);
		skelImg = cv::Mat(img.size(), CV_8UC1, cv::Scalar(0));
		cv::Mat temp;
		cv::Mat eroded;
		cv::Mat element = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3));
		bool done;
		do {
			cv::erode(img, eroded, element);
			cv::dilate(eroded, temp, element); // temp = open(img)
			cv::subtract(img, temp, temp);
			cv::bitwise_or(skelImg, temp, skelImg);
			eroded.copyTo(img);
			done = (cv::countNonZero(img) == 0);
		} while (!done);
	}
	// Skeletonization another method: ZS1984
	else {
		skeleton_thinning(img, skelImg);
	}

	// OUTPUT CV::POINTS
	//skelImg_out = skelImg.clone();

	cv::imshow("SkeletonWindonw", skelImg);
	cv::waitKey(0);

	if (is_solid)
		get_CVPoints_from_bwImg(skelImg, solid_skeleton_);
	else
		get_CVPoints_from_bwImg(skelImg, void_skeleton_);
}

void CSkeletonExtractor::get_CVPoints_from_bwImg(cv::Mat bwImg, std::vector<cv::Point>& pts)
{
	pts.clear();
	for (int c = 0; c < bwImg.cols; c++) {
		for (int r = 0; r < bwImg.rows; r++) {
			if (bwImg.at<uchar>(r, c) == 255) {
				pts.push_back(cv::Point(c, r));
			}
		}
	}
}

int CSkeletonExtractor::flood_fill_sampling(cv::Mat binaryGuide, std::vector<cv::Point>& samples, int radius)
{
	int num_criticalPts = 0;
	samples.clear();
	std::vector<std::pair<cv::Point, cv::Point>> sample_pair_list;
	cv::Mat bw = binaryGuide.clone();

	//cv::imshow("aa", binaryGuide);
	//cv::waitKey(0);

	// get bwPts
	std::vector<cv::Point> skeleton;
	get_CVPoints_from_bwImg(binaryGuide, skeleton);

	std::vector<cv::Point> joints, endpts;
	find_critical_pts(bw, skeleton, joints, endpts);
	std::vector<cv::Point> seeds;
	seeds.insert(seeds.end(), joints.begin(), joints.end());
	seeds.insert(seeds.end(), endpts.begin(), endpts.end());
	num_criticalPts = seeds.size();
	std::cout << "#joints: " << joints.size() << std::endl;
	std::cout << "#endpts: " << endpts.size() << std::endl;
	std::cout << "non zero count: " << cv::countNonZero(bw) << std::endl;

	// find critical points
	int niters = 0, max_iters = 100;
	do {
		// save seeds
		std::vector<cv::Point> cur_seeds = seeds;
		samples.insert(samples.end(), seeds.begin(), seeds.end());
		seeds.clear();

		// mask their neighborhoods
		cv::Mat mask = cv::Mat::ones(bw.size(), CV_8UC1);
		mask *= 255;
		for (int i = 0; i < cur_seeds.size(); i++)
		{
			std::vector<cv::Point> nns;
			square_element(cur_seeds[i], radius, nns);
			for (int j = 0; j < nns.size(); j++)
			{
				mask.at<uchar>(nns[j].y, nns[j].x) = 0;
			}
		}
		bw &= mask;
		//cv::imshow("edge", bw);
		//cv::waitKey(10);
		// find new seeds
		//int num_seed = 0;
		std::set<cv::Point, cvpt_compare> tmp_seed_set;
		for (int i = 0; i < cur_seeds.size(); i++)
		{
			int i_num_seed = 0;
			std::vector<cv::Point> bns;
			outer_skirt_element(cur_seeds[i], radius + 1, bns);
			for (int j = 0; j < bns.size(); j++)
			{
				if (bw.at<uchar>(bns[j].y, bns[j].x) == 255)
				{
					std::vector<cv::Point> nns;
					get_neighbors(bns[j], nns);
					int B = 0;
					for (int k = 0; k < nns.size() - 1; k++) {
						if (bw.at<uchar>(nns[k].y, nns[k].x) == 0 && bw.at<uchar>(nns[k + 1].y, nns[k + 1].x) == 255)
							B++;
					}
					if (nns.size() == 8) {
						if (bw.at<uchar>(nns[7].y, nns[7].x) == 0 && bw.at<uchar>(nns[0].y, nns[0].x) == 255)
							B++;
					}
					else if (bw.at<uchar>(nns[0].y, nns[0].x) == 255)
						B++;
					if (B != 2) {
						auto it = tmp_seed_set.insert(bns[j]);
						// if sucessfully inserted
						if (it.second) {
							sample_pair_list.push_back(std::make_pair(cur_seeds[i], bns[j]));
							//i_num_seed++;
						}
					}
				}
			}
		}
		seeds.insert(seeds.end(), tmp_seed_set.begin(), tmp_seed_set.end());

		//std::cout << "niter:" << niters << 
		//	", num seed found: " << seeds.size() <<
		//	//", num seed actual: " << joints.size() + endpts.size() <<
		//	", non zero count: " << cv::countNonZero(bw) << 
		//	//"\n "<< 
		//	std::endl;

		if ((niters++) >= max_iters) // must be at the end of the loop
			break;
	} while (cv::countNonZero(bw) > 0 && seeds.size() > 0);

	return num_criticalPts;
}

void CSkeletonExtractor::find_critical_pts(cv::Mat binaryGuide, std::vector<cv::Point> skeleton, 
	std::vector<cv::Point>& joints, std::vector<cv::Point>& endpts)
{
	joints.clear();
	endpts.clear();
	cv::Mat bwclone = binaryGuide.clone();
	//cv::threshold(bwclone, bwclone, 10, 1, CV_THRESH_BINARY);

	// Find joints
	for (int i = 0; i < skeleton.size(); i++)
	{
		int c = skeleton[i].x, r = skeleton[i].y;
		if (bwclone.at<uchar>(r, c) == 255) {
			std::vector<cv::Point> nns;
			get_neighbors(skeleton[i], nns);

			int B = 0; // substract itself
			for (int k = 0; k < nns.size() - 1; k++) {
				if (bwclone.at<uchar>(nns[k].y, nns[k].x) == 0 &&
					bwclone.at<uchar>(nns[k + 1].y, nns[k + 1].x) == 255)
				{
					B++;
				}
			}
			if (nns.size() == 8)
			{
				if (bwclone.at<uchar>(nns[7].y, nns[7].x) == 0 &&
					bwclone.at<uchar>(nns[0].y, nns[0].x) == 255)
				{
					B++;
				}
			}
			else if (bwclone.at<uchar>(nns[0].y, nns[0].x) == 255) {
				B++;
			}

			if (B > 2 || B == 0) {
				joints.push_back(cv::Point(c, r));
			}
			else if (B == 1)
			{
				endpts.push_back(cv::Point(c, r));
			}
		}
	}
}

void CSkeletonExtractor::get_neighbors(cv::Point Apt, std::vector<cv::Point>& nns)
{
	int width, height;
	img_data_.getBoundingDomain(width, height);

	int x = Apt.x;
	int y = Apt.y;
	if (x == 0)
	{
		if (y == 0) {
			nns.push_back(cv::Point(x, y + 1));
			nns.push_back(cv::Point(x + 1, y + 1));
			nns.push_back(cv::Point(x + 1, y));
		}
		else if (y == height - 1) {
			nns.push_back(cv::Point(x + 1, y));
			nns.push_back(cv::Point(x + 1, y - 1));
			nns.push_back(cv::Point(x, y - 1));
		}
		else {
			nns.push_back(cv::Point(x, y + 1));
			nns.push_back(cv::Point(x + 1, y + 1));
			nns.push_back(cv::Point(x + 1, y));
			nns.push_back(cv::Point(x + 1, y - 1));
			nns.push_back(cv::Point(x, y - 1));
		}
	}
	else if (x == width - 1)
	{
		if (y == 0) {
			nns.push_back(cv::Point(x - 1, y));
			nns.push_back(cv::Point(x - 1, y + 1));
			nns.push_back(cv::Point(x, y + 1));
		}
		else if (y == height - 1) {
			nns.push_back(cv::Point(x, y - 1));
			nns.push_back(cv::Point(x - 1, y - 1));
			nns.push_back(cv::Point(x - 1, y));

		}
		else {
			nns.push_back(cv::Point(x, y - 1));
			nns.push_back(cv::Point(x - 1, y - 1));
			nns.push_back(cv::Point(x - 1, y));
			nns.push_back(cv::Point(x - 1, y + 1));
			nns.push_back(cv::Point(x, y + 1));
		}
	}
	else {
		if (y == 0) {
			nns.push_back(cv::Point(x - 1, y));
			nns.push_back(cv::Point(x - 1, y + 1));
			nns.push_back(cv::Point(x, y + 1));
			nns.push_back(cv::Point(x + 1, y + 1));
			nns.push_back(cv::Point(x + 1, y));
		}
		else if (y == height - 1) {
			nns.push_back(cv::Point(x + 1, y));
			nns.push_back(cv::Point(x + 1, y - 1));
			nns.push_back(cv::Point(x, y - 1));
			nns.push_back(cv::Point(x - 1, y - 1));
			nns.push_back(cv::Point(x - 1, y));
		}
		else {
			nns.push_back(cv::Point(x - 1, y));
			nns.push_back(cv::Point(x - 1, y + 1));
			nns.push_back(cv::Point(x, y + 1));
			nns.push_back(cv::Point(x + 1, y + 1));
			nns.push_back(cv::Point(x + 1, y));
			nns.push_back(cv::Point(x + 1, y - 1));
			nns.push_back(cv::Point(x, y - 1));
			nns.push_back(cv::Point(x - 1, y - 1));
		}
	}
}


void CSkeletonExtractor::square_element(cv::Point Apt, int Rsize, std::vector<cv::Point>& nns)
{
	int width, height;
	img_data_.getBoundingDomain(width, height);

	int Nsize = 2 * Rsize + 1; // make sure this is an odd number
	for (int r = 0; r < Nsize; r++) {
		for (int c = 0; c < Nsize; c++) {
			int x = Apt.x + (c - Rsize);
			int y = Apt.y + (r - Rsize);
			if (x >= 0 && x < width && y >= 0 && y < height)
				nns.push_back(cv::Point(x, y));
		}
	}
}

void CSkeletonExtractor::outer_skirt_element(cv::Point Apt, int Rsize, std::vector<cv::Point>& nns)
{
	assert(Rsize > 1);

	int width, height;
	img_data_.getBoundingDomain(width, height);

	int Nsize = 2 * Rsize + 1; // make sure this is an odd number
	std::set<std::pair<int, int>> nn_set;
	int x, y;
	y = Apt.y - Rsize;
	if (y >= 0 && y < height) {
		for (int i = 0; i < Nsize; i++) {
			x = Apt.x + (i - Rsize);
			if (x >= 0 && x < width)
				nn_set.insert(std::make_pair(x, y));
		}
	}
	y = Apt.y + Rsize;
	if (y >= 0 && y < height) {
		for (int i = 0; i < Nsize; i++) {
			x = Apt.x + (i - Rsize);
			if (x >= 0 && x < width)
				nn_set.insert(std::make_pair(x, y));
		}
	}
	x = Apt.x - Rsize;
	if (x >= 0 && x < width) {
		for (int i = 0; i < Nsize; i++) {
			y = Apt.y + (i - Rsize);
			if (y >= 0 && y < height)
				nn_set.insert(std::make_pair(x, y));
		}
	}
	x = Apt.x + Rsize;
	if (x >= 0 && x < width) {
		for (int i = 0; i < Nsize; i++) {
			y = Apt.y + (i - Rsize);
			if (y >= 0 && y < height)
				nn_set.insert(std::make_pair(x, y));
		}
	}
	for (auto it_set = nn_set.begin(); it_set != nn_set.end(); ++it_set)
		nns.push_back(cv::Point(it_set->first, it_set->second));
}


void CSkeletonExtractor::skeleton_thinning(cv::Mat binaryInput, cv::Mat & binaryOutput)
{
	//cleanBoundary(binaryInput);
	// Zhang and Suen, 1984 (modified from https://github.com/arnaud-ramey/voronoi/blob/master/src/voronoi.h)
	int max_iters = 100;
	cv::Mat skel = binaryInput.clone();
	cv::threshold(skel, skel, 10, 1, cv::THRESH_BINARY);
	cv::Mat prev = cv::Mat::zeros(skel.size(), CV_8UC1);
	cv::Mat diff;
	int niters = 0;
	do {
		skeleton_thinning_ZS_iter(skel, 0);
		skeleton_thinning_ZS_iter(skel, 1);
		cv::absdiff(skel, prev, diff);
		skel.copyTo(prev);
		if ((niters++) >= max_iters) // must be at the end of the loop
			break;
	} while (cv::countNonZero(diff) > 0);
	// Output
	skel *= 255;
	binaryOutput = skel.clone();
}

void CSkeletonExtractor::skeleton_thinning_ZS_iter(cv::Mat & im, int iter)
{
	cv::Mat marker = cv::Mat::zeros(im.size(), CV_8UC1);
	for (int i = 1; i < im.rows - 1; i++)
	{
		for (int j = 1; j < im.cols - 1; j++)
		{
			uchar p2 = im.at<uchar>(i - 1, j);
			uchar p3 = im.at<uchar>(i - 1, j + 1);
			uchar p4 = im.at<uchar>(i, j + 1);
			uchar p5 = im.at<uchar>(i + 1, j + 1);
			uchar p6 = im.at<uchar>(i + 1, j);
			uchar p7 = im.at<uchar>(i + 1, j - 1);
			uchar p8 = im.at<uchar>(i, j - 1);
			uchar p9 = im.at<uchar>(i - 1, j - 1);

			int A = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
				(p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
				(p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
				(p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);
			int B = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
			int m1 = iter == 0 ? (p2 * p4 * p6) : (p2 * p4 * p8);
			int m2 = iter == 0 ? (p4 * p6 * p8) : (p2 * p6 * p8);

			if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
				marker.at<uchar>(i, j) = 1;
		}
	}
	im &= ~marker;
}

void CSkeletonExtractor::convert_skeleton_to_img(std::vector<cv::Point> skeleton, cv::Mat & img)
{
	int w, h;
	img_data_.getBoundingDomain(w, h);
	cv::Mat skel_img = cv::Mat(h, w, CV_8UC1);
	skel_img.setTo(0);
	
	for (int i = 0; i < skeleton.size(); i++)
	{
		cv::Point p = skeleton[i];
		skel_img.at<uchar>(p.y, p.x) = 255;
	}
	skel_img.copyTo(img);
}
